precision highp float;

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D iorTex;           // 4 (for metals)
uniform sampler2D kTex;             // 5 (for metals)
uniform sampler2D envMap;           // 6
in vec2 vTexCoord;

layout(location = 0) out vec4 gbuf_rad;
layout(location = 1) out vec4 gbuf_rng;

uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

uniform float minLengthScale;
uniform float maxLengthScale;
uniform float skyPower;
uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapRotation;
uniform float radianceClamp;
uniform float skipProbability;
uniform float shadowStrength;
uniform bool maxStepsIsMiss;

uniform vec3 surfaceDiffuseAlbedoRGB;

#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2
#define MAT_VOLUM  3

#define M_PI 3.1415926535897932384626433832795

//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

__DEFINES__

__SHADER__

///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit over specified segment
bool traceDistance(in vec3 start, in vec3 dir, float maxDist,
                   inout vec3 hit, inout int material)
{
    float minMarchDist = minLengthScale;

    const float HUGE_VAL = 1.0e20;
    float sdf = HUGE_VAL;
#ifdef HAS_SURFACE
    float sdf_surfa = abs(SDF_SURFACE(start));    sdf = min(sdf, sdf_surfa);
#endif
#ifdef HAS_METAL
    float sdf_metal = abs(SDF_METAL(start));      sdf = min(sdf, sdf_metal);
#endif
#ifdef HAS_DIELECTRIC
    float sdf_diele = abs(SDF_DIELECTRIC(start)); sdf = min(sdf, sdf_diele);
#endif
#ifdef HAS_VOLUME
    float sdf_volum = abs(SDF_VOLUME(start)); sdf = min(sdf, sdf_volum);
#endif

    float InitialSign = sign(sdf);
    float t = 0.0;
    int iters=0;
    for (int n=0; n<__MAX_MARCH_STEPS__; n++)
    {
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf;
        if (t>=maxDist) break;
        vec3 pW = start + t*dir;
        sdf = HUGE_VAL;
#ifdef HAS_SURFACE
        sdf_surfa = abs(SDF_SURFACE(pW));    if (sdf_surfa<minMarchDist) { material = MAT_SURFA; break; } sdf = min(sdf, sdf_surfa);
#endif
#ifdef HAS_METAL
        sdf_metal = abs(SDF_METAL(pW));      if (sdf_metal<minMarchDist) { material = MAT_METAL; break; } sdf = min(sdf, sdf_metal);
#endif
#ifdef HAS_DIELECTRIC
        sdf_diele = abs(SDF_DIELECTRIC(pW)); if (sdf_diele<minMarchDist) { material = MAT_DIELE; break; } sdf = min(sdf, sdf_diele); 
#endif
#ifdef HAS_VOLUME
        sdf_volum = abs(SDF_VOLUME(pW)); if (sdf_volum<minMarchDist) { material = MAT_VOLUM; break; } sdf = min(sdf, sdf_volum); 
#endif
        iters++;
    }
    hit = start + t*dir;
    if (t>=maxDist) return false;
    if (maxStepsIsMiss && iters>=__MAX_MARCH_STEPS__) return false;
    return true;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir,
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_INVAL;
    return traceDistance(start, dir, maxMarchDist, hit, material);
}

// (whether occluded along infinite ray)
bool Occluded(in vec3 start, in vec3 dir)
{
    float eps = 3.0*minLengthScale;
    vec3 delta = eps * dir;
    int material;
    vec3 hit;
    return traceRay(start+delta, dir, hit, material, maxLengthScale);
}

vec3 normal(in vec3 pW, int material)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = 2.0*minLengthScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 xyyp = pW+e.xyy; vec3 xyyn = pW-e.xyy;
    vec3 yxyp = pW+e.yxy; vec3 yxyn = pW-e.yxy;
    vec3 yyxp = pW+e.yyx; vec3 yyxn = pW-e.yyx;
    vec3 N;
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { N = vec3(   SDF_SURFACE(xyyp) -    SDF_SURFACE(xyyn),    SDF_SURFACE(yxyp) -    SDF_SURFACE(yxyn),    SDF_SURFACE(yyxp) -    SDF_SURFACE(yyxn)); return normalize(N); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { N = vec3(     SDF_METAL(xyyp) -      SDF_METAL(xyyn),      SDF_METAL(yxyp) -      SDF_METAL(yxyn),      SDF_METAL(yyxp) -      SDF_METAL(yyxn)); return normalize(N); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { N = vec3(SDF_DIELECTRIC(xyyp) - SDF_DIELECTRIC(xyyn), SDF_DIELECTRIC(yxyp) - SDF_DIELECTRIC(yxyn), SDF_DIELECTRIC(yyxp) - SDF_DIELECTRIC(yyxn)); return normalize(N); }
#endif
#ifdef HAS_VOLUME
    if (material==MAT_VOLUM) { N = vec3(    SDF_VOLUME(xyyp) -     SDF_VOLUME(xyyn),     SDF_VOLUME(yxyp) -     SDF_VOLUME(yxyn),     SDF_VOLUME(yyxp) -     SDF_VOLUME(yyxn)); return normalize(N); }
#endif
}

/////////////////////////////////////////////////////////////////////////
// Basis transforms
/////////////////////////////////////////////////////////////////////////

float cosTheta2(in vec3 nLocal) { return nLocal.z*nLocal.z; }
float cosTheta(in vec3 nLocal)  { return nLocal.z; }
float sinTheta2(in vec3 nLocal) { return 1.0 - cosTheta2(nLocal); }
float sinTheta(in vec3 nLocal)  { return sqrt(max(0.0, sinTheta2(nLocal))); }
float tanTheta2(in vec3 nLocal) { float ct2 = cosTheta2(nLocal); return max(0.0, 1.0 - ct2) / max(ct2, 1.0e-7); }
float tanTheta(in vec3 nLocal)  { return sqrt(max(0.0, tanTheta2(nLocal))); }

struct Basis
{
    vec3 nW; // aligned with the z-axis in local space
    vec3 tW;
    vec3 bW;
};

Basis makeBasis(in vec3 nW)
{
    Basis basis;
    basis.nW = nW;
    if (abs(nW.z) < abs(nW.x))
    {
        basis.tW.x =  nW.z;
        basis.tW.y =  0.0;
        basis.tW.z = -nW.x;
    }
    else
    {
        basis.tW.x =  0.0;
        basis.tW.y = nW.z;
        basis.tW.z = -nW.y;
    }
    basis.tW = normalize(basis.tW);
    basis.bW = cross(nW, basis.tW);
    return basis;
}

vec3 worldToLocal(in vec3 vWorld, in Basis basis)
{
    return vec3( dot(vWorld, basis.tW),
                 dot(vWorld, basis.bW),
                 dot(vWorld, basis.nW) );
}

vec3 localToWorld(in vec3 vLocal, in Basis basis)
{
    return basis.tW*vLocal.x + basis.bW*vLocal.y + basis.nW*vLocal.z;
}

vec3 xyzToRgb(vec3 XYZ)
{
    // (Assuming RGB in sRGB color space)
    vec3 RGB;
    RGB.r =  3.2404542*XYZ.x - 1.5371385*XYZ.y - 0.4985314*XYZ.z;
    RGB.g = -0.9692660*XYZ.x + 1.8760108*XYZ.y + 0.0415560*XYZ.z;
    RGB.b =  0.0556434*XYZ.x - 0.2040259*XYZ.y + 1.0572252*XYZ.z;
    return RGB;
}

vec3 rgbToXyz(in vec3 RGB)
{
    // (Assuming RGB in sRGB color space)
    vec3 XYZ;
    XYZ.x = 0.4124564*RGB.r + 0.3575761*RGB.g + 0.1804375*RGB.b;
    XYZ.y = 0.2126729*RGB.r + 0.7151522*RGB.g + 0.0721750*RGB.b;
    XYZ.z = 0.0193339*RGB.r + 0.1191920*RGB.g + 0.9503041*RGB.b;
    return XYZ;
}

vec3 xyz_to_spectrum(vec3 XYZ)
{
    // Given a color in XYZ tristimulus coordinates, return the coefficients in spectrum:
    //      L(l) = cx x(l) + cy y(l) + cz z(l)
    // (where x, y, z are the tristimulus CMFs) such that this spectrum reproduces the given XYZ tristimulus.
    // (NB, these coefficients differ from XYZ, because the XYZ color matching functions are not orthogonal)
    vec3 c;
    c.x =  3.38214566*XYZ.x - 2.58540997*XYZ.y - 0.40649004*XYZ.z;
    c.y = -2.58540997*XYZ.x + 3.20943158*XYZ.y + 0.22767094*XYZ.z;
    c.z = -0.40649004*XYZ.x + 0.22767094*XYZ.y + 0.70334476*XYZ.z;

    return c;
}

/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

/// GLSL floating point pseudorandom number generator, from
/// "Implementing a Photorealistic Rendering System using GLSL", Toshiya Hachisuka
/// http://arxiv.org/pdf/1505.06022.pdf
float rand(inout vec4 rnd)
{
    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);
    vec4 beta = floor(rnd/q);
    vec4 p = a*(rnd - beta*q) - beta*r;
    beta = (1.0 - sign(p))*0.5*m;
    rnd = p + beta;
    return fract(dot(rnd/m, vec4(1.0, -1.0, 1.0, -1.0)));
}

vec3 sampleHemisphere(inout vec4 rnd, inout float pdf)
{
    // Do cosine-weighted sampling of hemisphere
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}

vec3 SURFACE_DIFFUSE_REFL_RGB(in vec3 X, in vec3 nW, in vec3 woW)
{
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, woW);
    return reflRGB;
}

////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

vec3 environmentRadiance(in vec3 dir)
{
    float phi = atan(dir.x, dir.z) + M_PI + M_PI*envMapRotation/180.0;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = acos(dir.y);             // theta in [0, pi]
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 RGB;
    if (haveEnvMap)
    {
        RGB = texture(envMap, vec2(u,v)).rgb;
        RGB *= skyPower;
    }
    else
    {
        RGB = skyPower * vec3(1.0);
    }
    return RGB;
}

////////////////////////////////////////////////////////////////////////////////
// Ambient occlusion integrator
////////////////////////////////////////////////////////////////////////////////

void constructPrimaryRay(in vec2 pixel, inout vec4 rnd,
                         inout vec3 primaryStart, inout vec3 primaryDir)
{
    // Compute world ray direction for given (possibly jittered) fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
    float fh = tan(0.5*radians(camFovy)); // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    primaryDir = normalize(camDir + s);
    if (camAperture<=0.0)
    {
        primaryStart = camPos;
        return;
    }
    vec3 focalPlaneHit = camPos + camFocalDistance*primaryDir/dot(primaryDir, camDir);
    float lensRadial = camAperture * sqrt(rand(rnd));
    float theta = 2.0*M_PI * rand(rnd);
    vec3 lensPos = camPos + lensRadial*(-camX*cos(theta) + camY*sin(theta));
    primaryStart = lensPos;
    primaryDir = normalize(focalPlaneHit - lensPos);
}

void main()
{
    INIT();

    vec4 rnd = texture(RngData, vTexCoord);
    if (rand(rnd) < skipProbability)
    {
        vec4 oldL = texture(Radiance, vTexCoord);
        float oldN = oldL.w;
        float newN = oldN;
        vec3 newL = oldL.rgb;

        gbuf_rad = vec4(newL, newN);
        gbuf_rng = rnd;
        return;
    }

    // Jitter over pixel
    vec2 pixel = gl_FragCoord.xy;
    pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    vec3 primaryStart, primaryDir;
    constructPrimaryRay(pixel, rnd, primaryStart, primaryDir);

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int hitMaterial;
    bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxLengthScale);

    vec3 RGB;
    if (hit)
    {
        // Compute normal at hit point
        vec3 nW = normal(pW, hitMaterial);
        Basis basis = makeBasis(nW);

        // Construct a uniformly sampled AO ray
        float hemispherePdf;
        vec3 wiL = sampleHemisphere(rnd, hemispherePdf);
        vec3 wiW = localToWorld(wiL, basis);

        // Compute diffuse albedo
        vec3 diffuseAlbedo = SURFACE_DIFFUSE_REFL_RGB(pW, basis.nW, woW);
        
        // Set incident radiance to according to whether the AO ray hit anything or missed.
        if (!Occluded(pW, wiW)) RGB = diffuseAlbedo;
        else                    RGB = diffuseAlbedo * abs(1.0 - shadowStrength);
    }
    else
    {
        if (envMapVisible) RGB = environmentRadiance(primaryDir);
        else               RGB = vec3(0.0);
    }
    
    vec3 XYZ = rgbToXyz(RGB);
    
    // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + XYZ) / newN;

    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}
