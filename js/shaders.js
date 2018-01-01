var Shaders = {

'ao-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D RadianceBlocks;   // 4
uniform sampler2D iorTex;           // 5 (for metals)
uniform sampler2D kTex;             // 6 (for metals)
uniform sampler2D envMap;           // 7
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
uniform bool jitter;

uniform float metalRoughness;
uniform vec3 metalSpecAlbedoXYZ;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoXYZ;
uniform vec3 surfaceDiffuseAlbedoXYZ;
uniform vec3 surfaceSpecAlbedoXYZ;
uniform float surfaceRoughness;
uniform float surfaceIor;

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
    vec3 C = xyzToRgb(surfaceDiffuseAlbedoXYZ);
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(C, X, nW, woW);
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
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

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
`,

'ao-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'firsthit-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D RadianceBlocks;   // 4
uniform sampler2D iorTex;           // 5 (for metals)
uniform sampler2D kTex;             // 6 (for metals)
uniform sampler2D envMap;           // 7
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
uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;

uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapRotation;
uniform float radianceClamp;
uniform float skipProbability;
uniform float shadowStrength;
uniform bool maxStepsIsMiss;
uniform bool jitter;

uniform vec3 surfaceDiffuseAlbedoXYZ;
uniform vec3 surfaceSpecAlbedoXYZ;
uniform float surfaceRoughness;
uniform float surfaceIor;

#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5

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
bool traceDistance(in vec3 start, in vec3 dir, float maxDist, inout vec3 hit)
{
    float minMarch = minLengthScale;
    const float HUGE_VAL = 1.0e20;
    float sdf = abs(SDF_SURFACE(start));
    float InitialSign = sign(sdf);
    float t = InitialSign * sdf; // (always take the first step along the ray direction)
    if (t>=maxDist) return false;
    for (int n=0; n<__MAX_MARCH_STEPS__; n++)
    {
        vec3 pW = start + t*dir;
        float sdf_surfa = abs(SDF_SURFACE(pW));
        if (sdf_surfa<minMarch)
        {
            hit = start + t*dir;
            return true;
        }
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf_surfa;
        if (t>=maxDist) return false;
    }
    return !maxStepsIsMiss;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir, inout vec3 hit, float maxMarchDist)
{
    return traceDistance(start, dir, maxMarchDist, hit);
}

vec3 normal(in vec3 pW)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = 2.0*minLengthScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 xyyp = pW+e.xyy; vec3 xyyn = pW-e.xyy;
    vec3 yxyp = pW+e.yxy; vec3 yxyn = pW-e.yxy;
    vec3 yyxp = pW+e.yyx; vec3 yyxn = pW-e.yyx;
    vec3 N = vec3(SDF_SURFACE(xyyp)-SDF_SURFACE(xyyn), SDF_SURFACE(yxyp)-SDF_SURFACE(yxyn), SDF_SURFACE(yyxp)- SDF_SURFACE(yyxn));
    return normalize(N);
}

// MC-estimates the amount of light transmitted along an infinite ray
float Visibility(in vec3 pW, in vec3 rayDir)
{
    vec3 pW_next;
    bool hit = traceRay(pW, rayDir, pW_next, maxLengthScale);
    if (!hit) return 1.0; // if no hit, return visibility 1.0 (transparent)
    return 0.0;
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
    Basis b;
    b.nW = nW;
    if (abs(nW.z) < abs(nW.x)) { b.tW.x =  nW.z; b.tW.y =  0.0; b.tW.z = -nW.x; }
    else                       { b.tW.x =  0.0;  b.tW.y = nW.z; b.tW.z = -nW.y; }
    b.tW = normalize(b.tW);
    b.bW = cross(nW, b.tW);
    return b;
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


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

/// Uniformly sample a sphere
vec3 sampleSphere(inout vec4 rnd, inout float pdf)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphere(inout vec4 rnd, inout float pdf)
{
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}

// pdf for cosine-weighted sampling of hemisphere
float pdfHemisphere(in vec3 wiL)
{
    return wiL.z / M_PI;
}

float powerHeuristic(const float a, const float b)
{
    return a/(a + b);
}

vec3 samplePhaseFunction(in vec3 rayDir, inout vec4 rnd)
{
    float g = 0.0; // @todo: parameter
    float costheta;
    if (abs(g)<1.0e-3)
        costheta = 1.0 - 2.0*rand(rnd);
    else
        costheta = 1.0/(2.0*g) * (1.0 + g*g - ((1.0-g*g)*(1.0-g+2.0*g*rand(rnd))));
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0*M_PI*rand(rnd);
    Basis basis = makeBasis(rayDir);
    return costheta*rayDir + sintheta*(cos(phi)*basis.tW + sin(phi)*basis.bW);
}

float phaseFunction(float mu)
{
    float g = 0.0; // @todo: parameter
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}


/////////////////////////////////////////////////////////////////////////
// Beckmann Microfacet formulae
/////////////////////////////////////////////////////////////////////////

// m = the microfacet normal (in the local space where z = the macrosurface normal)
float microfacetEval(in vec3 m, in float roughness)
{
    float t2 = tanTheta2(m);
    float c2 = cosTheta2(m);
    float roughnessSqr = roughness*roughness;
    float epsilon = 1.0e-9;
    float exponent = t2 / max(roughnessSqr, epsilon);
    float D = exp(-exponent) / max(M_PI*roughnessSqr*c2*c2, epsilon);
    return D;
}

// m = the microfacet normal (in the local space where z = the macrosurface normal)
vec3 microfacetSample(inout vec4 rnd, in float roughness)
{
    float phiM = (2.0 * M_PI) * rand(rnd);
    float cosPhiM = cos(phiM);
    float sinPhiM = sin(phiM);
    float epsilon = 1.0e-9;
    float tanThetaMSqr = -roughness*roughness * log(max(epsilon, rand(rnd)));
    float cosThetaM = 1.0 / sqrt(1.0 + tanThetaMSqr);
    float sinThetaM = sqrt(max(0.0, 1.0 - cosThetaM*cosThetaM));
    return normalize(vec3(sinThetaM*cosPhiM, sinThetaM*sinPhiM, cosThetaM));
}

float microfacetPDF(in vec3 m, in float roughness)
{
    return microfacetEval(m, roughness) * abs(cosTheta(m));
}

// Shadow-masking function
// Approximation from Walter et al (v = arbitrary direction, m = microfacet normal)
float smithG1(in vec3 vLocal, in vec3 mLocal, float roughness)
{
    float tanThetaAbs = abs(tanTheta(vLocal));
    float epsilon = 1.0e-6;
    if (tanThetaAbs < epsilon) return 1.0; // perpendicular incidence -- no shadowing/masking
    if (dot(vLocal, mLocal) * vLocal.z <= 0.0) return 0.0; // Back side is not visible from the front side, and the reverse.
    float a = 1.0 / (max(roughness, epsilon) * tanThetaAbs); // Rational approximation to the shadowing/masking function (Walter et al)  (<0.35% rel. error)
    if (a >= 1.6) return 1.0;
    float aSqr = a*a;
    return (3.535*a + 2.181*aSqr) / (1.0 + 2.276*a + 2.577*aSqr);
}

float smithG2(in vec3 woL, in vec3 wiL, in vec3 mLocal, float roughness)
{
    return smithG1(woL, mLocal, roughness) * smithG1(wiL, mLocal, roughness);
}


/////////////////////////////////////////////////////////////////////////
// BRDF functions
/////////////////////////////////////////////////////////////////////////

vec3 SURFACE_DIFFUSE_REFL_RGB(in vec3 X, in vec3 nW, in vec3 woW)
{
    vec3 C = xyzToRgb(surfaceDiffuseAlbedoXYZ);
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(C, X, nW, woW);
    return reflRGB;
}

vec3 SURFACE_SPEC_REFL_RGB(in vec3 X, in vec3 nW, in vec3 woW)
{
    vec3 C = xyzToRgb(surfaceSpecAlbedoXYZ);
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(C, X, nW, woW);
    return reflRGB;
}

// Fast path for non-transmissive surface
float fresnelDielectricReflectanceFast(in float cosi, in float ior)
{
    float ei = 1.0;
    float et = max(1.0, ior);
    float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));
    float cost = sqrt(max(0.0, 1.0 - sint*sint));
    float cosip = abs(cosi);
    float rParallel      = (et*cosip - ei*cost) / (et*cosip + ei*cost);
    float rPerpendicular = (ei*cosip - et*cost) / (ei*cosip + et*cost);
    return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}

vec3 evaluateSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 diffuseAlbedo = SURFACE_DIFFUSE_REFL_RGB(X, basis.nW, woW);
    vec3    specAlbedo = SURFACE_SPEC_REFL_RGB(X, basis.nW, woW);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float specWeight = length(specAlbedo);
    float diffWeight = length(diffuseAlbedo);
    float weightSum = max(specWeight + diffWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    vec3 f = Fr * specAlbedo * specProb * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    f += diffuseAlbedo * (1.0-specProb)/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 diffuseAlbedo = SURFACE_DIFFUSE_REFL_RGB(X, basis.nW, woW);
    vec3    specAlbedo = SURFACE_SPEC_REFL_RGB(X, basis.nW, woW);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float diffusePdf = pdfHemisphere(wiL);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    float specWeight = length(specAlbedo);
    float diffWeight = length(diffuseAlbedo);
    float weightSum = max(specWeight + diffWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

vec3 sampleSurface(in vec3 X, in Basis basis, in vec3 woL, inout vec3 wiL, inout float pdfOut, inout vec4 rnd)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 diffuseAlbedo = SURFACE_DIFFUSE_REFL_RGB(X, basis.nW, woW);
    vec3    specAlbedo = SURFACE_SPEC_REFL_RGB(X, basis.nW, woW);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float specWeight = length(specAlbedo);
    float diffWeight = length(diffuseAlbedo);
    float weightSum = max(specWeight + diffWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        wiL = sampleHemisphere(rnd, pdfOut);
        pdfOut *= (1.0-specProb);
        return diffuseAlbedo*diffWeight/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
        if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
        vec3 f = Fr * specAlbedo * specWeight * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = specProb * microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
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

Basis sunBasis;

vec3 sampleSunDir(inout vec4 rnd)
{
    float costhetamax = cos(sunAngularSize * M_PI/180.0);
    float costheta = 1.0 - rand(rnd)*(1.0-costhetamax);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    return localToWorld(vec3(x, y, z), sunBasis);
}

vec3 sunRadiance(in vec3 dir)
{
    if (dot(dir, sunDir) < cos(sunAngularSize*M_PI/180.0)) return vec3(0.0);
    return sunPower * sunColor;
}

vec3 sampleLightAtSurface(Basis basis, inout vec4 rnd, inout vec3 wiL, inout vec3 wiW, inout float lightPdf)
{
    // Light sampling (choose either sun or sky)
    vec3 wiW_sun = sampleSunDir(rnd);
    vec3 wiL_sun = worldToLocal(wiW_sun, basis);
    float sunPdf = 1.0; // convenient to decouple total sun power from its angular size
    float sunWeight = sunPower * max(0.0, wiL_sun.z);

    float skyPdf;
    vec3 wiL_sky = sampleHemisphere(rnd, skyPdf);
    float skyWeight = skyPower * max(0.0, wiL_sky.z);

    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    bool chooseSun = (rand(rnd) <= sunProb);
    if (chooseSun)
    {
        lightPdf = sunPdf * sunProb;
        wiL = wiL_sun;
        wiW = wiW_sun;
        return sunRadiance(wiW_sun);
    }
    lightPdf = skyPdf * max(PDF_EPSILON, 1.0-sunProb);
    wiL = wiL_sky;
    wiW = localToWorld(wiL_sky, basis);
    return environmentRadiance(wiW);
}

// Estimate direct radiance at the given surface vertex
vec3 directSurfaceLighting(in vec3 pW, Basis basis, in vec3 woW, inout vec4 rnd, inout float lightPdf)
{
    vec3 wiL, wiW; // direction of sampled direct light (*towards* the light)
    vec3 Li = sampleLightAtSurface(basis, rnd, wiL, wiW, lightPdf);
    vec3 dPw = 3.0 * minLengthScale * basis.nW;
    float V = Visibility(pW+dPw, wiW);
    Li *= abs(1.0 - shadowStrength*(1.0-V));

    // Apply MIS weight with the BSDF pdf for the sampled direction
    vec3 woL = worldToLocal(woW, basis);
    float bsdfPdf = pdfSurface(pW, basis, woL, wiL);
    if (bsdfPdf<PDF_EPSILON) return vec3(0.0);
    vec3 f = evaluateSurface(pW, basis, woL, wiL);
    float misWeight = powerHeuristic(lightPdf, bsdfPdf);
    vec3 fOverPdf = min(vec3(radianceClamp), f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li * abs(dot(wiW, basis.nW)) * misWeight;
}


////////////////////////////////////////////////////////////////////////////////
// First hit integrator
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
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    vec3 primaryStart, primaryDir;
    constructPrimaryRay(pixel, rnd, primaryStart, primaryDir);

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    bool hit = traceRay(primaryStart, primaryDir, pW, maxLengthScale);

    // @todo: features to support realistic 'fractal architecture' renderings:
    //        - optionally dim according to visibility with less dimming at greater hit distance
    //        - optional single volume scattering assuming homogeneous medium
    //               for fast god rays and distance fog.

    vec3 RGB = vec3(0.0, 0.0, 0.0);
    if (hit)
    {
        vec3 nW = normal(pW);
        Basis basis = makeBasis(nW);

        // Sample BSDF for the bounce direction
        vec3 woL = worldToLocal(woW, basis);
        vec3 wiL;
        float bsdfPdf;
        vec3 f = sampleSurface(pW, basis, woL, wiL, bsdfPdf, rnd);

        vec3 wiW = localToWorld(wiL, basis);
        vec3 rayDir = wiW; // bounce ray direction

        // Update path throughput
        vec3 fOverPdf = min(vec3(radianceClamp), f/max(PDF_EPSILON, bsdfPdf));
        vec3 throughput = fOverPdf * abs(dot(wiW, nW));

        // Add term from direct light sampling:
        float lightPdf;
        RGB += throughput * directSurfaceLighting(pW, basis, woW, rnd, lightPdf);

        // Add lighting term from BRDF sampling (optionally)
        if (__MAX_BOUNCES__>0)
        {
            float eps = 3.0*minLengthScale;
            pW += nW * sign(dot(wiW, nW)) * eps;
            float V = Visibility(pW, rayDir);
            vec3 Li = environmentRadiance(rayDir);
            Li += sunRadiance(rayDir);
            Li *= abs(1.0 - shadowStrength*(1.0-V));
            RGB += throughput * Li * powerHeuristic(bsdfPdf, lightPdf);
        }
    }
    else
    {
        if (envMapVisible)      RGB += environmentRadiance(primaryDir);
        if (sunVisibleDirectly) RGB += sunRadiance(primaryDir);
    }

    vec3 XYZ = rgbToXyz(RGB);

    // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.xyz + XYZ) / newN;

    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}
`,

'firsthit-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'normals-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D RadianceBlocks;   // 4
uniform sampler2D iorTex;           // 5 (for metals)
uniform sampler2D kTex;             // 6 (for metals)
uniform sampler2D envMap;           // 7
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
uniform bool jitter;

uniform float metalRoughness;
uniform vec3 metalSpecAlbedoXYZ;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoXYZ;
uniform vec3 surfaceDiffuseAlbedoXYZ;
uniform vec3 surfaceSpecAlbedoXYZ;
uniform float surfaceRoughness;
uniform float surfaceIor;

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

vec3 rgbToXyz(in vec3 RGB)
{
    // (Assuming RGB in sRGB color space)
    vec3 XYZ;
    XYZ.x = 0.4124564*RGB.r + 0.3575761*RGB.g + 0.1804375*RGB.b;
    XYZ.y = 0.2126729*RGB.r + 0.7151522*RGB.g + 0.0721750*RGB.b;
    XYZ.z = 0.0193339*RGB.r + 0.1191920*RGB.g + 0.9503041*RGB.b;
    return XYZ;
}

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

////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

vec3 environmentRadianceXYZ(in vec3 dir)
{
    float phi = atan(dir.x, dir.z) + M_PI + M_PI*envMapRotation/180.0;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = acos(dir.y);           // theta in [0, pi]
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 XYZ;
    if (haveEnvMap)
    {
        vec3 RGB = texture(envMap, vec2(u,v)).rgb;
        RGB *= skyPower;
        XYZ = rgbToXyz(RGB);
    }
    else
    {
        XYZ = skyPower * vec3(1.0);
    }
    return XYZ;
}

////////////////////////////////////////////////////////////////////////////////
// Normals integrator
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
    vec2 pixel = gl_FragCoord.xy;
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    vec3 primaryStart, primaryDir;
    constructPrimaryRay(pixel, rnd, primaryStart, primaryDir);

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int hitMaterial;
    bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxLengthScale);
    vec3 colorXYZ;
    if (hit)
    {
        // Compute normal at hit point
        vec3 nW = normal(pW, hitMaterial);
        colorXYZ = rgbToXyz(0.5*(nW+vec3(1.0)));
    }
    else
    {
        if (envMapVisible) colorXYZ = environmentRadianceXYZ(primaryDir);
        else               colorXYZ = vec3(0.0);
    }


   // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;

    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}
`,

'normals-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'pathtracer-fragment-shader-tmp': `#version 300 es
precision highp float;

vec3 test(in vec3 pW, int material);

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D RadianceBlocks;   // 4
uniform sampler2D iorTex;           // 4
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

uniform float lengthScale;
uniform float minLengthScale;
uniform float maxLengthScale;

uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;

uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapRotation;
uniform float skyPower;

uniform float radianceClamp;
uniform float skipProbability;
uniform float shadowStrength;
uniform bool maxStepsIsMiss;
uniform bool jitter;
uniform int wavelengthSamples;

uniform float metalRoughness;
uniform vec3 metalSpecAlbedoXYZ;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoXYZ;
uniform vec3 surfaceDiffuseAlbedoXYZ;
uniform vec3 surfaceSpecAlbedoXYZ;
uniform float surfaceRoughness;
uniform float surfaceIor;

uniform float volumeExtinction;
uniform vec3 volumeScatteringColorRGB;
uniform vec3 volumeAbsorptionColorRGB;
uniform float volumeEmission;
uniform vec3 volumeEmissionColorRGB;
uniform float volumeAnisotropy;

#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2
#define MAT_VOLUM  3

#define M_PI 3.1415926535897932384626433832795

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

//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

__DEFINES__

__SHADER__

__IOR_FUNC__


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

Basis sunBasis;

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

vec3 rgbToSpectrum(vec3 RGB)
{
    vec3 SPC;
    SPC.x =  0.66656354*RGB.x + 0.20530493*RGB.y + 0.06824062*RGB.z;
    SPC.y = -1.94519051*RGB.x + 2.95689669*RGB.y + 0.21337152*RGB.z;
    SPC.z = -1.05576383*RGB.x + 0.92124612*RGB.y + 1.15412969*RGB.z;
    return SPC;
}


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

/// Uniformly sample a sphere
vec3 sampleSphere(inout vec4 rnd, inout float pdf)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphere(inout vec4 rnd, inout float pdf)
{
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}

// pdf for cosine-weighted sampling of hemisphere
float pdfHemisphere(in vec3 wiL)
{
    return wiL.z / M_PI;
}

float powerHeuristic(const float a, const float b)
{
    return a/(a + b);
}

vec3 samplePhaseFunction(in vec3 rayDir, inout vec4 rnd)
{
    float g = volumeAnisotropy;
    float costheta;
    if (abs(g)<1.0e-3)
        costheta = 1.0 - 2.0*rand(rnd);
    else
        costheta = 1.0/(2.0*g) * (1.0 + g*g - ((1.0-g*g)*(1.0-g+2.0*g*rand(rnd))));
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0*M_PI*rand(rnd);
    Basis basis = makeBasis(rayDir);
    return costheta*rayDir + sintheta*(cos(phi)*basis.tW + sin(phi)*basis.bW);
}

float phaseFunction(float mu)
{
    float g = volumeAnisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}


/////////////////////////////////////////////////////////////////////////
// Beckmann Microfacet formulae
/////////////////////////////////////////////////////////////////////////

// m = the microfacet normal (in the local space where z = the macrosurface normal)
float microfacetEval(in vec3 m, in float roughness)
{
    float t2 = tanTheta2(m);
    float c2 = cosTheta2(m);
    float roughnessSqr = roughness*roughness;
    float epsilon = 1.0e-9;
    float exponent = t2 / max(roughnessSqr, epsilon);
    float D = exp(-exponent) / max(M_PI*roughnessSqr*c2*c2, epsilon);
    return D;
}

// m = the microfacet normal (in the local space where z = the macrosurface normal)
vec3 microfacetSample(inout vec4 rnd, in float roughness)
{
    float phiM = (2.0 * M_PI) * rand(rnd);
    float cosPhiM = cos(phiM);
    float sinPhiM = sin(phiM);
    float epsilon = 1.0e-9;
    float tanThetaMSqr = -roughness*roughness * log(max(epsilon, rand(rnd)));
    float cosThetaM = 1.0 / sqrt(1.0 + tanThetaMSqr);
    float sinThetaM = sqrt(max(0.0, 1.0 - cosThetaM*cosThetaM));
    return normalize(vec3(sinThetaM*cosPhiM, sinThetaM*sinPhiM, cosThetaM));
}

float microfacetPDF(in vec3 m, in float roughness)
{
    return microfacetEval(m, roughness) * abs(cosTheta(m));
}

// Shadow-masking function
// Approximation from Walter et al (v = arbitrary direction, m = microfacet normal)
float smithG1(in vec3 vLocal, in vec3 mLocal, float roughness)
{
    float tanThetaAbs = abs(tanTheta(vLocal));
    float epsilon = 1.0e-6;
    if (tanThetaAbs < epsilon) return 1.0; // perpendicular incidence -- no shadowing/masking
    if (dot(vLocal, mLocal) * vLocal.z <= 0.0) return 0.0; // Back side is not visible from the front side, and the reverse.
    float a = 1.0 / (max(roughness, epsilon) * tanThetaAbs); // Rational approximation to the shadowing/masking function (Walter et al)  (<0.35% rel. error)
    if (a >= 1.6) return 1.0;
    float aSqr = a*a;
    return (3.535*a + 2.181*aSqr) / (1.0 + 2.276*a + 2.577*aSqr);
}

float smithG2(in vec3 woL, in vec3 wiL, in vec3 mLocal, float roughness)
{
    return smithG1(woL, mLocal, roughness) * smithG1(wiL, mLocal, roughness);
}


/////////////////////////////////////////////////////////////////////////
// BRDF functions
/////////////////////////////////////////////////////////////////////////

#ifdef HAS_GEOMETRY

// ****************************        Dielectric        ****************************

/// Compute Fresnel reflectance at a dielectric interface (which has an "interior" and an "exterior").
/// Here cosi is the cosine to the (interior-to-exterior) normal of the incident ray
/// direction wi, (where we use the PBRT convention that the light
/// travels in the opposite direction to wi, i.e. wi is the direction
/// *from* which the light arrives at the surface point).
/// The Boolean "entering" indicates whether the ray is entering or exiting the dielectric.
float fresnelDielectricReflectance(in float cosi, in float iorInternal, in float iorExternal)
{
    bool entering = cosi > 0.0;
    float ei, et;
    if (entering)
    {
        ei = iorExternal; // Incident from external medium, if entering
        et = iorInternal; // Transmitted to internal medium, if entering
    }
    else
    {
        ei = iorInternal; // Incident from internal medium, if exiting
        et = iorExternal; // Transmitted to external medium, if exiting
    }

    // Compute sint from Snell's law:
    float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));

    // Handle total internal reflection
    if (sint >= 1.0) return 1.0;

    float cost = sqrt(max(0.0, 1.0 - sint*sint));
    float cosip = abs(cosi);
    float rParallel      = (et*cosip - ei*cost) / (et*cosip + ei*cost);
    float rPerpendicular = (ei*cosip - et*cost) / (ei*cosip + et*cost);
    return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}

// Given the direction (wt) of a beam transmitted through a plane dielectric interface
// with the given normal (n) (with n pointing away from the transmitted half-space, i.e. wt.n<0),
// and the ratio eta = et/ei of the incident IOR (ei) and transmitted IOR (et),
// compute the direction of the incident beam (wi).
// Returns false if no such beam exists, due to total internal reflection.
bool refraction(in vec3 n, in float eta, in vec3 wt, inout vec3 wi)
{
    vec3 x = normalize(cross(cross(n, wt), n));
    float sint = dot(wt, x);
    float sini = eta * sint; // Snell's law
    float sini2 = sini*sini;
    if (sini2 >= 1.0) return false; // Total internal reflection
    float cosi2 = (1.0 - sini*sini);
    float cosi = max(0.0, sqrt(cosi2));
    wi = -cosi*n + sini*x;
    return true;
}

float DIELECTRIC_SPEC_REFL_EVAL(in vec3 X, in vec3 woL, in Basis basis, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 C = xyzToRgb(dieleSpecAlbedoXYZ); // pass dieleSpecAlbedoRGB instead
    vec3 reflRGB = DIELECTRIC_SPECULAR_REFLECTANCE(C, X, basis.nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

float evaluateDielectric( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    float dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 h;
    float eta; // IOR ratio, et/ei
    if (reflected)
    { // Compute reflection half-vector
        h = normalize(wiL + woL);
    }
    else
    { // Compute refraction half-vector
        bool entering = cosTheta(wiL) > 0.0;
        eta = entering ? ior : 1.0/ior;
        h = normalize(wiL + eta*woL);
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f;
    if (reflected)
    {
        f = dielectricAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    }
    else
    {
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        float dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        f = (1.0 - Fr) * G * D * abs(im) * dwh_dwo / max(abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    }
    return f;
}

float pdfDielectric( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    float dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 h;
    float dwh_dwo;
    float pdf;
    if (reflected)
    {
        h = normalize(wiL + woL);
        dwh_dwo = 1.0 / max(4.0*dot(woL, h), DENOM_TOLERANCE);
        pdf = Fr;
    }
    else
    { // Compute reflection half-vector
        bool entering = cosTheta(wiL) > 0.0;
        float eta = entering ? ior : 1.0/max(ior, DENOM_TOLERANCE);
        vec3 h = normalize(wiL + eta*woL);
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        pdf = 1.0 - Fr;
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    pdf *= microfacetPDF(h, roughness);
    return abs(pdf * dwh_dwo);
}

float sampleDielectric( in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 XYZ,
                        inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_DIELE(wavelength_nm);
    float dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    float microPDF = microfacetPDF(m, roughness);
    float reflectProb = Fr;
    if (rand(rnd) < reflectProb)
    {
        // Compute specularly reflected ray direction
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute incident direction by reflecting woL about m
        if (wiL.z<0.0) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float f = Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = microPDF * reflectProb * dwh_dwo; // Return total BRDF and corresponding pdf
        return f;
    }
    else
    {
        // Note, for transmission:
        //  woL.m < 0 means the transmitted light is entering the dielectric *from* the (vacuum) exterior of the microfacet
        //  woL.m > 0 means the transmitted light is exiting the dielectric *to* the (vacuum) exterior of the microfacet
        bool entering = dot(woL, m) < 0.0;
        float eta;
        vec3 ni; // normal pointing into incident halfspace
        if (entering) // Entering, incident halfspace is outside dielectric
        {
            eta = ior; // = et/ei
            ni = m;
        }
        else // Exiting, incident halfspace is inside dielectric
        {
            eta = 1.0/ior; // = et/ei
            ni = -m;
        }
        // Compute incident direction corresponding to known transmitted direction
        if ( !refraction(ni, eta, woL, wiL) ) return 0.0; // total internal reflection occurred
        wiL = -wiL; // As refract() computes the incident beam direction, and wiL is defined to be opposite to that.
        float cosi = dot(wiL, m);
        float Frm = fresnelDielectricReflectance(cosi, ior, 1.0);
        float Trm = 1.0 - Frm;  // Fresnel transmittance
        // Evaluate microfacet distribution for the sampled half direction
        vec3 wh = m; // refraction half-vector = m
        float D = microfacetEval(wh, roughness);
        float G = smithG2(woL, wiL, wh, roughness); // Shadow-masking function
        float dwh_dwo; // Jacobian of the half-direction mapping
        float im = dot(wiL, m);
        {
            float om = dot(woL, m);
            float sqrtDenom = im + eta*om;
            dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        }
        float f = abs(im) * dwh_dwo * Trm * G * D / max(abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        pdfOut = (1.0-reflectProb) * microPDF * abs(dwh_dwo);
        return f;
    }
}

// ****************************        Metal        ****************************

/// cosi is the cosine to the (outward) normal of the incident ray direction wi,
/// ior is the index of refraction of the metal, and k its absorption coefficient
float fresnelMetalReflectance(in float cosi, in float ior, in float k)
{
    float cosip = abs(cosi);
    float cosi2 = cosip * cosip;
    float iikk  = ior*ior + k*k;
    float iikkc = iikk * cosi2;
    float twoEtaCosi = 2.0*ior*cosip;
    float Rparl2 = (iikkc - twoEtaCosi + 1.0  ) / (iikkc + twoEtaCosi + 1.0  );
    float Rperp2 = (iikk  - twoEtaCosi + cosi2) / (iikk  + twoEtaCosi + cosi2);
    return 0.5*(Rparl2 + Rperp2);
}

float METAL_SPEC_REFL_EVAL(in vec3 X, in vec3 woL, in Basis basis, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 C = xyzToRgb(metalSpecAlbedoXYZ); // pass metalSpecAlbedoRGB instead
    vec3 reflRGB = METAL_SPECULAR_REFLECTANCE(C, X, basis.nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

float evaluateMetal( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(woL.z, ior, k);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float specAlbedo = METAL_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    return f;
}

float pdfMetal( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float pdf = microfacetPDF(h, roughness) * dwh_dwo;
    return pdf;
}

float sampleMetal( in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 XYZ,
                   inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(woL.z, ior, k);
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
    if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
    float D = microfacetEval(m, roughness);
    float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
    float specAlbedo = METAL_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    float dwh_dwo; // Jacobian of the half-direction mapping
    dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
    pdfOut = microfacetPDF(m, roughness) * dwh_dwo;
    return f;
}


// ****************************        Surface        ****************************

float SURFACE_DIFFUSE_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 woW, in vec3 XYZ)
{
    vec3 C = xyzToRgb(surfaceDiffuseAlbedoXYZ); // pass surfaceDiffuseAlbedoRGB instead
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(C, X, nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

float SURFACE_SPEC_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 woW, in vec3 XYZ)
{
    vec3 C = xyzToRgb(surfaceSpecAlbedoXYZ); // pass surfaceSpecAlbedoRGB instead
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(C, X, nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

// Fast path for non-transmissive surface
float fresnelDielectricReflectanceFast(in float cosi, in float ior)
{
    float ei = 1.0;
    float et = max(1.0, ior);
    float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));
    float cost = sqrt(max(0.0, 1.0 - sint*sint));
    float cosip = abs(cosi);
    float rParallel      = (et*cosip - ei*cost) / (et*cosip + ei*cost);
    float rPerpendicular = (ei*cosip - et*cost) / (ei*cosip + et*cost);
    return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}

float evaluateSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    float diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, XYZ);
    float    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, XYZ);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    f += diffuseAlbedo/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    float diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, XYZ);
    float    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, XYZ);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float diffusePdf = pdfHemisphere(wiL);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    float sum = max(Fr*specAlbedo + diffuseAlbedo, DENOM_TOLERANCE);
    float specProb = Fr*specAlbedo/sum;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

float sampleSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 XYZ,
                    inout vec3 wiL, inout float pdfOut, inout vec4 rnd)
{
    vec3 woW = localToWorld(woL, basis);
    float diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, XYZ);
    float    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, XYZ);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float sum = max(specAlbedo + diffuseAlbedo, DENOM_TOLERANCE);
    float specProb = specAlbedo/sum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        wiL = sampleHemisphere(rnd, pdfOut);
        pdfOut *= (1.0-specProb);
        return diffuseAlbedo/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
        if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
        float f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = specProb * microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
}

// ****************************        BSDF common interface        ****************************

float evaluateBsdf( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ,
                    inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    evaluateSurface(X, basis, woL, wiL,                XYZ); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      evaluateMetal(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return evaluateDielectric(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
}

float sampleBsdf( in vec3 X, in Basis basis, in vec3 woL, in int material, in float wavelength_nm, in vec3 XYZ,
                  inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    sampleSurface(X, basis, woL,                XYZ, wiL, pdfOut, rnd); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      sampleMetal(X, basis, woL, wavelength_nm, XYZ, wiL, pdfOut, rnd); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return sampleDielectric(X, basis, woL, wavelength_nm, XYZ, wiL, pdfOut, rnd); }
#endif
}

float pdfBsdf( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    pdfSurface(X, basis, woL, wiL,                XYZ); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      pdfMetal(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return pdfDielectric(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
}

#endif // HAS_GEOMETRY

#ifdef HAS_VOLUME

float VOLUME_ALBEDO_EVAL(in vec3 X, in vec3 XYZ)
{
    vec3 sigma_s = VOLUME_SCATTERING_COLOR(volumeScatteringColorRGB, X);
    vec3 sigma_a = VOLUME_ABSORPTION_COLOR(volumeAbsorptionColorRGB, X);
    vec3 albedo = sigma_s / max(sigma_s + sigma_a, 1.0e-8);
    vec3 albedoSPC = rgbToSpectrum(reflRGB);
    return max(dot(albedoSPC, XYZ), 0.0);
}

float VOLUME_EXTINCTION_EVAL(in vec3 X, in vec3 XYZ)
{
    vec3 sigma_s = VOLUME_SCATTERING_COLOR(volumeScatteringColorRGB, X);
    vec3 sigma_a = VOLUME_ABSORPTION_COLOR(volumeAbsorptionColorRGB, X);
    float extinction = volumeExtinction * VOLUME_EXTINCTION(volumeExtinction, X);
    vec3 sigma_t = extinction * (sigma_s + sigma_a) / lengthScale;
    vec3 sigmaSPC = rgbToSpectrum(sigma_t);
    return max(0.0, dot(sigmaSPC, XYZ));
}

float VOLUME_EXTINCTION_MAX_EVAL()
{
    return VOLUME_EXTINCTION_MAX(volumeExtinction);
}

float VOLUME_EMISSION_EVAL(in vec3 X, in vec3 XYZ)
{
    vec3 emission = VOLUME_EMISSION(volumeEmissionColorRGB * volumeEmission, X);
    vec3 emissionSPC = rgbToSpectrum(emission);
    return max(0.0, dot(emissionSPC, XYZ));
}

float VOLUME_ANISOTROPY_EVAL(in vec3 X)
{
    return VOLUME_ANISOTROPY(volumeAnisotropy, X);
}

#endif


///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit over specified segment
bool traceDistance(in vec3 start, in vec3 dir, float maxDist,
                   inout vec3 hit, inout int material)
{
    float minMarch = minLengthScale;
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
    float sdf_volum = abs(SDF_VOLUME(start));     sdf = min(sdf, sdf_volum);
#endif

    float InitialSign = sign(sdf);
    float t = InitialSign * sdf; // (always take the first step along the ray direction)
    if (t>=maxDist) return false;
    for (int n=0; n<__MAX_MARCH_STEPS__; n++)
    {
        vec3 pW = start + t*dir;
        sdf = HUGE_VAL;
#ifdef HAS_SURFACE
        sdf_surfa = abs(SDF_SURFACE(pW));    if (sdf_surfa<minMarch) { material = MAT_SURFA; hit = start + t*dir; return true; } sdf = min(sdf, sdf_surfa);
#endif
#ifdef HAS_METAL
        sdf_metal = abs(SDF_METAL(pW));      if (sdf_metal<minMarch) { material = MAT_METAL; hit = start + t*dir; return true; } sdf = min(sdf, sdf_metal);
#endif
#ifdef HAS_DIELECTRIC
        sdf_diele = abs(SDF_DIELECTRIC(pW)); if (sdf_diele<minMarch) { material = MAT_DIELE; hit = start + t*dir; return true; } sdf = min(sdf, sdf_diele);
#endif
#ifdef HAS_VOLUME
        sdf_volum = abs(SDF_VOLUME(pW));     if (sdf_volum<minMarch) { material = MAT_VOLUM; hit = start + t*dir; return true; } sdf = min(sdf, sdf_volum);
#endif
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf;
        if (t>=maxDist) return false;
    }
    return !maxStepsIsMiss;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir,
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_INVAL;
    return traceDistance(start, dir, maxMarchDist, hit, material);
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

// MC-estimates the amount of light transmitted along an infinite ray
float Visibility(in vec3 pW, in vec3 rayDir, in vec3 XYZ, inout vec4 rnd, bool inVolume)
{
    vec3 pW_next;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
#ifndef HAS_VOLUME
    return float(!hit);
#else
    // Deal with all volume interactions until the ray either hits a surface or escapes
    const int maxVolHits = 16;
    int volHits = 0;
    float sigma_t_max = VOLUME_EXTINCTION_MAX_EVAL();
    float inv_sigma_t_max = 1.0/sigma_t_max;
    float eps = 3.0*minLengthScale;
    while ((inVolume || hitMaterial==MAT_VOLUM) && volHits<maxVolHits)
    {
        // If ray lies in the volume interior
        if (inVolume)
        {
            // Do Woodcock tracking to check for a possible volume scattering event along the in-volume segment pW -> pW_next
            // (except if the current ray lies inside a dielectric, in which case there is no volumetric scattering, i.e. the dielectric displaces the volume)
            vec3 pScatter;
            float sigma_t;
            bool interacted = false;
            float segmentLength = length(pW_next - pW);
            float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max;
            int steps = 0;
            while (distanceMarched<segmentLength-eps && steps<__MAX_VOLUME_STEPS__)
            {
                pScatter = pW + distanceMarched*rayDir;
                sigma_t = VOLUME_EXTINCTION_EVAL(pScatter, XYZ);
                interacted = bool(rand(rnd) < sigma_t * inv_sigma_t_max);
                if (interacted) break;
                distanceMarched += -log(rand(rnd)) * inv_sigma_t_max;
                steps++;
            }

            if (interacted) return 0.0; // on absorption/scattering, return zero visibility
            else
            {
                pW = pW_next;
                if (hitMaterial==MAT_VOLUM) // if hitMaterial is volume boundary, displace pW into volume exterior
                {
                    vec3 nW = normal(pW, MAT_VOLUM);
                    pW += nW * eps;
                    inVolume = false;
                }
                else // if hitMaterial is a surface, return zero visibility
                    return 0.0;
            }
        }

        // else if ray passed through the volume exterior and hit the volume boundary
        else // !inVolume
        {
            pW = pW_next; // advance to the volume boundary
            vec3 nW = normal(pW, MAT_VOLUM); // displace pW into volume interior
            pW -= nW * eps;
            inVolume = true;
        }

        // Re-trace from new pW
        hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
        volHits++;
    }
    return float(!hit || volHits==maxVolHits); // if no hit or scattering, return visibility 1.0 (transparent)
#endif // HAS_VOLUME
}


////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dir)
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

float environmentRadiance(in vec3 dir, in vec3 XYZ)
{
    vec3 RGB_sky = environmentRadianceRGB(dir);
    vec3 XYZ_spec = rgbToSpectrum(RGB_sky); // convert to radiance at the given wavelength
    return dot(XYZ, XYZ_spec);
}

vec3 sampleSunDir(inout vec4 rnd)
{
    float costhetamax = cos(sunAngularSize * M_PI/180.0);
    float costheta = 1.0 - rand(rnd)*(1.0-costhetamax);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    return localToWorld(vec3(x, y, z), sunBasis);
}

float sunRadiance(in vec3 dir, in vec3 XYZ)
{
    if (dot(dir, sunDir) < cos(sunAngularSize*M_PI/180.0)) return 0.0;
    vec3 RGB_sun = sunPower * sunColor;
    vec3 XYZ_spec = rgbToXyz(RGB_sun); // convert to radiance at the given wavelength
    return dot(XYZ, XYZ_spec);
}

#ifdef HAS_GEOMETRY

float sampleLightAtSurface(Basis basis, in vec3 XYZ, inout vec4 rnd, inout vec3 wiL, inout vec3 wiW, inout float lightPdf)
{
    // Light sampling (choose either sun or sky)
    vec3 wiW_sun = sampleSunDir(rnd);
    vec3 wiL_sun = worldToLocal(wiW_sun, basis);
    float sunPdf = 1.0; // convenient to decouple total sun power from its angular size
    float sunWeight = sunPower * max(0.0, wiL_sun.z);

    float skyPdf;
    vec3 wiL_sky = sampleHemisphere(rnd, skyPdf);
    float skyWeight = skyPower * max(0.0, wiL_sky.z);

    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    bool chooseSun = (rand(rnd) <= sunProb);
    if (chooseSun)
    {
        lightPdf = sunPdf * sunProb;
        wiL = wiL_sun;
        wiW = wiW_sun;
        return sunRadiance(wiW_sun, XYZ);
    }
    lightPdf = skyPdf * max(PDF_EPSILON, 1.0-sunProb);
    wiL = wiL_sky;
    wiW = localToWorld(wiL_sky, basis);
    return environmentRadiance(wiW, XYZ);
}

// Estimate direct radiance at the given surface vertex
float directSurfaceLighting(in vec3 pW, Basis basis, in vec3 woW, in int material,
                            float wavelength_nm, in vec3 XYZ, inout vec4 rnd, bool inVolume, inout float lightPdf)
{
    vec3 wiL, wiW; // direction of sampled direct light (*towards* the light)
    float Li = sampleLightAtSurface(basis, XYZ, rnd, wiL, wiW, lightPdf);
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    float V = Visibility(pW+dPw, wiW, XYZ, rnd, inVolume);
    Li *= abs(1.0 - shadowStrength*(1.0-V));

    // Apply MIS weight with the BSDF pdf for the sampled direction
    vec3 woL = worldToLocal(woW, basis);
    float bsdfPdf = pdfBsdf(pW, basis, woL, wiL, material, wavelength_nm, XYZ);
    if (bsdfPdf<PDF_EPSILON) return 0.0;
    float f = evaluateBsdf(pW, basis, woL, wiL, material, wavelength_nm, XYZ, rnd);
    float misWeight = powerHeuristic(lightPdf, bsdfPdf);
    float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li * abs(dot(wiW, basis.nW)) * misWeight;
}

#endif // HAS_GEOMETRY

#ifdef HAS_VOLUME

float sampleLightInVolume(in vec3 XYZ, inout vec4 rnd, inout vec3 wiW, inout float lightPdf)
{
    // Light sampling (choose either sun or sky)
    vec3 wiW_sun = sampleSunDir(rnd);
    float sunPdf = 1.0; // convenient to decouple total sun power from its angular size
    float sunWeight = sunPower * sunPdf;

    float skyPdf;
    vec3 wiW_sky = sampleSphere(rnd, skyPdf);
    float skyWeight = skyPower * skyPdf;

    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    bool chooseSun = (rand(rnd) <= sunProb);
    if (chooseSun)
    {
        lightPdf = sunPdf * sunProb;
        wiW = wiW_sun;
        return sunRadiance(wiW_sun, XYZ);
    }
    lightPdf = skyPdf * max(PDF_EPSILON, 1.0-sunProb);
    wiW = wiW_sky;
    return environmentRadiance(wiW, XYZ);
}

// Estimate direct radiance at the given volumetric vertex
float directVolumeLighting(in vec3 pW, in vec3 woW, in vec3 XYZ, inout vec4 rnd, bool inVolume)
{
    vec3 wiW; // direction of sampled direct light (*towards* the light)
    float lightPdf;
    float Li = sampleLightInVolume(XYZ, rnd, wiW, lightPdf);
    float V = Visibility(pW, wiW, XYZ, rnd, inVolume);
    Li *= abs(1.0 - shadowStrength*(1.0-V));

    float f = phaseFunction(dot(woW, -wiW));
    float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li;
}

#endif // HAS_VOLUME

////////////////////////////////////////////////////////////////////////////////
// Pathtracing integrator
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
    float theta = 2.0*M_PI*rand(rnd);
    vec3 lensPos = camPos + lensRadial*(-camX*cos(theta) + camY*sin(theta));
    primaryStart = lensPos;
    primaryDir = normalize(focalPlaneHit - lensPos);
}

void pathtrace(vec2 pixel, vec4 rnd) // the current pixel
{
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

    // Sample photon wavelength via the inverse CDF of the emission spectrum.
    // We limit the sampled wavelengths to between about 420nm and 700nm to avoid color mismatch at low sample count
    float xi = 0.1+0.8*(0.5+floor(float(wavelengthSamples)*rand(rnd))) / float(wavelengthSamples);
    float w = texture(ICDF, vec2(xi, 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;

    // Convert wavelength to XYZ tristimulus
    vec3 XYZ = texture(WavelengthToXYZ, vec2(w, 0.5)).rgb;
    vec3 RGB = clamp(xyzToRgb(XYZ), 0.0, 1.0);
    
    // Jitter over pixel
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    // Compute world ray direction for this fragment
    vec3 primaryStart, primaryDir;
    constructPrimaryRay(pixel, rnd, primaryStart, primaryDir);

    // Perform pathtrace to estimate the primary ray radiance, L
    float L = 0.0;
    float throughput = 1.0;
    float misWeight = 1.0; // For MIS book-keeping
    vec3 pW = primaryStart;
    vec3 rayDir = primaryDir; // (opposite to light direction)
    bool inVolume = false;
#ifdef HAS_VOLUME
    float sigma_t_max = VOLUME_EXTINCTION_MAX_EVAL();
    float inv_sigma_t_max = 1.0/sigma_t_max;
    inVolume = SDF_VOLUME(primaryStart) < 0.0;
#endif
    float eps = 3.0*minLengthScale;
#ifdef HAS_DIELECTRIC
    bool inDielectric = SDF_DIELECTRIC(primaryStart) < 0.0;
#endif

    for (int vertex=0; vertex<=__MAX_BOUNCES__; ++vertex)
    {
        // Raycast along current propagation direction rayDir, from current vertex pW to pW_next
        vec3 pW_next;
        int hitMaterial;
        bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
        float rayLength = maxLengthScale;
        if (hit)
        {
            rayLength = length(pW_next - pW);
        }

#ifdef HAS_VOLUME
        // Deal with all volume interactions until the ray either hits a surface or escapes
        int volHits = 0;
        bool absorb = false;
        while ((inVolume || hitMaterial==MAT_VOLUM) && volHits<(1+__MAX_SCATTERS__))
        {
            // If ray lies in the volume interior
            if (inVolume)
            {
                float segmentLength = length(pW_next - pW);

                // Do Woodcock tracking to check for a possible volume absorption/scattering event along the in-volume segment pW -> pW_next
                // (except if the current ray lies inside a dielectric, in which case there is no volumetric scattering, i.e. the dielectric displaces the volume)
                vec3 pScatter;
                float sigma_t;
                bool interacted = false;
#ifdef HAS_DIELECTRIC
                if (!inDielectric)
#endif
                {
                    float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max;
                    int steps = 0;
                    while (distanceMarched<segmentLength-eps && steps<__MAX_VOLUME_STEPS__)
                    {
                        pScatter = pW + distanceMarched*rayDir;
                        sigma_t = VOLUME_EXTINCTION_EVAL(pScatter, XYZ);
                        interacted = bool(rand(rnd) < sigma_t * inv_sigma_t_max);
                        if (interacted) break;
                        distanceMarched += -log(rand(rnd)) * inv_sigma_t_max;
                        steps++;
                    }
                }

                // Absorption/scattering event
                if (interacted)
                {
                    vec3 wiW = samplePhaseFunction(rayDir, rnd); // sample scattered dir (NB, PF is the angle PDF so the MC PDF denom. cancels it in the L estimator)
                    float emission = VOLUME_EMISSION_EVAL(pScatter, XYZ);
                    float albedo = VOLUME_ALBEDO_EVAL(pScatter, XYZ);
                    absorb = (rand(rnd) > albedo);
                    if (absorb)
                    {
                        L += throughput * emission; // add volume emission term
                        break;
                    }

                    L += throughput * directVolumeLighting(pScatter, -rayDir, XYZ, rnd, inVolume); // add contribution due to direct lighting at vertex
                    throughput *= albedo; // update throughput due to scattering
                    if (throughput < THROUGHPUT_EPSILON) break;
                    pW = pScatter; // continue to next vertex
                    rayDir = wiW;
                }

                // No scattering
                else
                {
                    pW = pW_next;
                    if (hitMaterial==MAT_VOLUM) // if hitMaterial is volume boundary, displace pW into volume exterior
                    {
                        vec3 nW = normal(pW, MAT_VOLUM);
                        pW += nW * eps;
                        inVolume = false;
                    }
                    else        // if hitMaterial is a surface, break and handle it (note, still inside volume)
                        break;  // (or if !hit, that's a numerical glitch as we should be inside a closed volume boundary)
                }
            }

            // else if ray passed through the volume exterior and hit the volume boundary
            else // !inVolume
            {
                pW = pW_next; // advance to the volume boundary
                vec3 nW = normal(pW, MAT_VOLUM); // displace pW into volume interior
                pW -= nW * eps;
                inVolume = true;
            }

            // Re-trace from new pW
            hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
            volHits++;
        }

        // If ray was absorbed by the volume, terminate
        if (absorb) break;
        if (volHits == __MAX_SCATTERS__) hit = false;
#endif // HAS_VOLUME

        // If this ray didn't scatter and missed all geometry, add environment light term and terminate path
        if (!hit)
        {
            float Li = 0.0;
            bool sunVisible = true; // @todo: parameter
            if (!(vertex==0 && !envMapVisible))      Li += environmentRadiance(rayDir, XYZ);
            if (!(vertex==0 && !sunVisibleDirectly)) Li += sunRadiance(rayDir, XYZ);
            L += throughput * Li * misWeight;
            break;
        }

#ifndef HAS_GEOMETRY
        break;
#else
        // If the current ray lies inside a dielectric, apply Beer's law for absorption
#ifdef HAS_DIELECTRIC
        if (inDielectric)
        {
            throughput *= exp(-rayLength*dot(dieleAbsorptionRGB, RGB));
        }
#endif

        // This ray didn't scatter but hit some geometry, so deal with the surface interaction.
        // First, compute normal at current surface vertex
        pW = pW_next;
        vec3 nW = normal(pW, hitMaterial);
        Basis basis = makeBasis(nW);

        // Sample BSDF for the next bounce direction
        vec3 woW = -rayDir;
        vec3 woL = worldToLocal(woW, basis);
        vec3 wiL;
        float bsdfPdf;
        float f = sampleBsdf(pW, basis, woL, hitMaterial, wavelength_nm, XYZ, wiL, bsdfPdf, rnd);
        vec3 wiW = localToWorld(wiL, basis);
        rayDir = wiW; // Update ray direction

        // Detect dielectric transmission
#ifdef HAS_DIELECTRIC
        if (hitMaterial==MAT_DIELE && dot(woW, nW)*dot(wiW, nW) < 0.0)
        {
            inDielectric = !inDielectric;
        }
#endif

        // Update path throughput
        float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, bsdfPdf));
        throughput *= fOverPdf * abs(dot(wiW, nW));
        if (throughput < THROUGHPUT_EPSILON) break;

        // Add direct lighting term at current surface vertex
        float lightPdf = 0.0;
#ifdef HAS_DIELECTRIC
        if (!inDielectric)
#endif
        {
            L += throughput * directSurfaceLighting(pW, basis, woW, hitMaterial, wavelength_nm, XYZ, rnd, inVolume, lightPdf);
        }

        // Prepare for tracing the bounce ray
        pW += nW * sign(dot(wiW, nW)) * eps;           // perturb vertex into half-space of scattered ray
        misWeight = powerHeuristic(bsdfPdf, lightPdf); // compute MIS weight for bounce ray
#endif // HAS_GEOMETRY
    }

    // Compute tristimulus contribution from estimated radiance
    vec3 colorXYZ = XYZ * L;

    // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;

    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}

void main()
{
    INIT();
    vec4 rnd = texture(RngData, vTexCoord);
    pathtrace(gl_FragCoord.xy, rnd);
}
`,

'pathtracer-fragment-shader': `#version 300 es
precision highp float;

vec3 test(in vec3 pW, int material);

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D RadianceBlocks;   // 4
uniform sampler2D iorTex;           // 4
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

uniform float lengthScale;
uniform float minLengthScale;
uniform float maxLengthScale;

uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;

uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapRotation;
uniform float skyPower;

uniform float radianceClamp;
uniform float skipProbability;
uniform float shadowStrength;
uniform bool maxStepsIsMiss;
uniform bool jitter;
uniform int wavelengthSamples;

uniform float metalRoughness;
uniform vec3 metalSpecAlbedoXYZ;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoXYZ;
uniform vec3 surfaceDiffuseAlbedoXYZ;
uniform vec3 surfaceSpecAlbedoXYZ;
uniform float surfaceRoughness;
uniform float surfaceIor;

uniform float volumeExtinction;
uniform vec3 volumeScatteringColorRGB;
uniform vec3 volumeAbsorptionColorRGB;
uniform float volumeEmission;
uniform vec3 volumeEmissionColorRGB;
uniform float volumeAnisotropy;

#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2
#define MAT_VOLUM  3

#define M_PI 3.1415926535897932384626433832795

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

//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

__DEFINES__

__SHADER__

__IOR_FUNC__


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

Basis sunBasis;

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

vec3 rgbToSpectrum(vec3 RGB)
{
    vec3 SPC;
    SPC.x =  1.48660973*RGB.x - 1.22162347*RGB.y  + 0.06635964*RGB.z;
    SPC.y = -0.67364598*RGB.x + 2.48197587*RGB.y  - 0.03286052*RGB.z;
    SPC.z = -0.1875686*RGB.x  + 0.17986238 *RGB.y + 1.08568973*RGB.z;
    return SPC;
}


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

/// Uniformly sample a sphere
vec3 sampleSphere(inout vec4 rnd, inout float pdf)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphere(inout vec4 rnd, inout float pdf)
{
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}

// pdf for cosine-weighted sampling of hemisphere
float pdfHemisphere(in vec3 wiL)
{
    return wiL.z / M_PI;
}

float powerHeuristic(const float a, const float b)
{
    return a/(a + b);
}

vec3 samplePhaseFunction(in vec3 rayDir, inout vec4 rnd)
{
    float g = volumeAnisotropy;
    float costheta;
    if (abs(g)<1.0e-3)
        costheta = 1.0 - 2.0*rand(rnd);
    else
        costheta = 1.0/(2.0*g) * (1.0 + g*g - ((1.0-g*g)*(1.0-g+2.0*g*rand(rnd))));
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0*M_PI*rand(rnd);
    Basis basis = makeBasis(rayDir);
    return costheta*rayDir + sintheta*(cos(phi)*basis.tW + sin(phi)*basis.bW);
}

float phaseFunction(float mu)
{
    float g = volumeAnisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}


/////////////////////////////////////////////////////////////////////////
// Beckmann Microfacet formulae
/////////////////////////////////////////////////////////////////////////

// m = the microfacet normal (in the local space where z = the macrosurface normal)
float microfacetEval(in vec3 m, in float roughness)
{
    float t2 = tanTheta2(m);
    float c2 = cosTheta2(m);
    float roughnessSqr = roughness*roughness;
    float epsilon = 1.0e-9;
    float exponent = t2 / max(roughnessSqr, epsilon);
    float D = exp(-exponent) / max(M_PI*roughnessSqr*c2*c2, epsilon);
    return D;
}

// m = the microfacet normal (in the local space where z = the macrosurface normal)
vec3 microfacetSample(inout vec4 rnd, in float roughness)
{
    float phiM = (2.0 * M_PI) * rand(rnd);
    float cosPhiM = cos(phiM);
    float sinPhiM = sin(phiM);
    float epsilon = 1.0e-9;
    float tanThetaMSqr = -roughness*roughness * log(max(epsilon, rand(rnd)));
    float cosThetaM = 1.0 / sqrt(1.0 + tanThetaMSqr);
    float sinThetaM = sqrt(max(0.0, 1.0 - cosThetaM*cosThetaM));
    return normalize(vec3(sinThetaM*cosPhiM, sinThetaM*sinPhiM, cosThetaM));
}

float microfacetPDF(in vec3 m, in float roughness)
{
    return microfacetEval(m, roughness) * abs(cosTheta(m));
}

// Shadow-masking function
// Approximation from Walter et al (v = arbitrary direction, m = microfacet normal)
float smithG1(in vec3 vLocal, in vec3 mLocal, float roughness)
{
    float tanThetaAbs = abs(tanTheta(vLocal));
    float epsilon = 1.0e-6;
    if (tanThetaAbs < epsilon) return 1.0; // perpendicular incidence -- no shadowing/masking
    if (dot(vLocal, mLocal) * vLocal.z <= 0.0) return 0.0; // Back side is not visible from the front side, and the reverse.
    float a = 1.0 / (max(roughness, epsilon) * tanThetaAbs); // Rational approximation to the shadowing/masking function (Walter et al)  (<0.35% rel. error)
    if (a >= 1.6) return 1.0;
    float aSqr = a*a;
    return (3.535*a + 2.181*aSqr) / (1.0 + 2.276*a + 2.577*aSqr);
}

float smithG2(in vec3 woL, in vec3 wiL, in vec3 mLocal, float roughness)
{
    return smithG1(woL, mLocal, roughness) * smithG1(wiL, mLocal, roughness);
}


/////////////////////////////////////////////////////////////////////////
// BRDF functions
/////////////////////////////////////////////////////////////////////////

#ifdef HAS_GEOMETRY

// ****************************        Dielectric        ****************************

/// Compute Fresnel reflectance at a dielectric interface (which has an "interior" and an "exterior").
/// Here cosi is the cosine to the (interior-to-exterior) normal of the incident ray
/// direction wi, (where we use the PBRT convention that the light
/// travels in the opposite direction to wi, i.e. wi is the direction
/// *from* which the light arrives at the surface point).
/// The Boolean "entering" indicates whether the ray is entering or exiting the dielectric.
float fresnelDielectricReflectance(in float cosi, in float iorInternal, in float iorExternal)
{
    bool entering = cosi > 0.0;
    float ei, et;
    if (entering)
    {
        ei = iorExternal; // Incident from external medium, if entering
        et = iorInternal; // Transmitted to internal medium, if entering
    }
    else
    {
        ei = iorInternal; // Incident from internal medium, if exiting
        et = iorExternal; // Transmitted to external medium, if exiting
    }

    // Compute sint from Snell's law:
    float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));

    // Handle total internal reflection
    if (sint >= 1.0) return 1.0;

    float cost = sqrt(max(0.0, 1.0 - sint*sint));
    float cosip = abs(cosi);
    float rParallel      = (et*cosip - ei*cost) / (et*cosip + ei*cost);
    float rPerpendicular = (ei*cosip - et*cost) / (ei*cosip + et*cost);
    return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}

// Given the direction (wt) of a beam transmitted through a plane dielectric interface
// with the given normal (n) (with n pointing away from the transmitted half-space, i.e. wt.n<0),
// and the ratio eta = et/ei of the incident IOR (ei) and transmitted IOR (et),
// compute the direction of the incident beam (wi).
// Returns false if no such beam exists, due to total internal reflection.
bool refraction(in vec3 n, in float eta, in vec3 wt, inout vec3 wi)
{
    vec3 x = normalize(cross(cross(n, wt), n));
    float sint = dot(wt, x);
    float sini = eta * sint; // Snell's law
    float sini2 = sini*sini;
    if (sini2 >= 1.0) return false; // Total internal reflection
    float cosi2 = (1.0 - sini*sini);
    float cosi = max(0.0, sqrt(cosi2));
    wi = -cosi*n + sini*x;
    return true;
}

float DIELECTRIC_SPEC_REFL_EVAL(in vec3 X, in vec3 woL, in Basis basis, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 C = xyzToRgb(dieleSpecAlbedoXYZ); // pass dieleSpecAlbedoRGB instead
    vec3 reflRGB = DIELECTRIC_SPECULAR_REFLECTANCE(C, X, basis.nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

float evaluateDielectric( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    float dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 h;
    float eta; // IOR ratio, et/ei
    if (reflected)
    { // Compute reflection half-vector
        h = normalize(wiL + woL);
    }
    else
    { // Compute refraction half-vector
        bool entering = cosTheta(wiL) > 0.0;
        eta = entering ? ior : 1.0/ior;
        h = normalize(wiL + eta*woL);
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f;
    if (reflected)
    {
        f = dielectricAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    }
    else
    {
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        float dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        f = (1.0 - Fr) * G * D * abs(im) * dwh_dwo / max(abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    }
    return f;
}

float pdfDielectric( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    float dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 h;
    float dwh_dwo;
    float pdf;
    if (reflected)
    {
        h = normalize(wiL + woL);
        dwh_dwo = 1.0 / max(4.0*dot(woL, h), DENOM_TOLERANCE);
        pdf = Fr;
    }
    else
    { // Compute reflection half-vector
        bool entering = cosTheta(wiL) > 0.0;
        float eta = entering ? ior : 1.0/max(ior, DENOM_TOLERANCE);
        vec3 h = normalize(wiL + eta*woL);
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        pdf = 1.0 - Fr;
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    pdf *= microfacetPDF(h, roughness);
    return abs(pdf * dwh_dwo);
}

float sampleDielectric( in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 XYZ,
                        inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_DIELE(wavelength_nm);
    float dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    float microPDF = microfacetPDF(m, roughness);
    float reflectProb = Fr;
    if (rand(rnd) < reflectProb)
    {
        // Compute specularly reflected ray direction
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute incident direction by reflecting woL about m
        if (wiL.z<0.0) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float f = Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = microPDF * reflectProb * dwh_dwo; // Return total BRDF and corresponding pdf
        return f;
    }
    else
    {
        // Note, for transmission:
        //  woL.m < 0 means the transmitted light is entering the dielectric *from* the (vacuum) exterior of the microfacet
        //  woL.m > 0 means the transmitted light is exiting the dielectric *to* the (vacuum) exterior of the microfacet
        bool entering = dot(woL, m) < 0.0;
        float eta;
        vec3 ni; // normal pointing into incident halfspace
        if (entering) // Entering, incident halfspace is outside dielectric
        {
            eta = ior; // = et/ei
            ni = m;
        }
        else // Exiting, incident halfspace is inside dielectric
        {
            eta = 1.0/ior; // = et/ei
            ni = -m;
        }
        // Compute incident direction corresponding to known transmitted direction
        if ( !refraction(ni, eta, woL, wiL) ) return 0.0; // total internal reflection occurred
        wiL = -wiL; // As refract() computes the incident beam direction, and wiL is defined to be opposite to that.
        float cosi = dot(wiL, m);
        float Frm = fresnelDielectricReflectance(cosi, ior, 1.0);
        float Trm = 1.0 - Frm;  // Fresnel transmittance
        // Evaluate microfacet distribution for the sampled half direction
        vec3 wh = m; // refraction half-vector = m
        float D = microfacetEval(wh, roughness);
        float G = smithG2(woL, wiL, wh, roughness); // Shadow-masking function
        float dwh_dwo; // Jacobian of the half-direction mapping
        float im = dot(wiL, m);
        {
            float om = dot(woL, m);
            float sqrtDenom = im + eta*om;
            dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        }
        float f = abs(im) * dwh_dwo * Trm * G * D / max(abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        pdfOut = (1.0-reflectProb) * microPDF * abs(dwh_dwo);
        return f;
    }
}

// ****************************        Metal        ****************************

/// cosi is the cosine to the (outward) normal of the incident ray direction wi,
/// ior is the index of refraction of the metal, and k its absorption coefficient
float fresnelMetalReflectance(in float cosi, in float ior, in float k)
{
    float cosip = abs(cosi);
    float cosi2 = cosip * cosip;
    float iikk  = ior*ior + k*k;
    float iikkc = iikk * cosi2;
    float twoEtaCosi = 2.0*ior*cosip;
    float Rparl2 = (iikkc - twoEtaCosi + 1.0  ) / (iikkc + twoEtaCosi + 1.0  );
    float Rperp2 = (iikk  - twoEtaCosi + cosi2) / (iikk  + twoEtaCosi + cosi2);
    return 0.5*(Rparl2 + Rperp2);
}

float METAL_SPEC_REFL_EVAL(in vec3 X, in vec3 woL, in Basis basis, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 C = xyzToRgb(metalSpecAlbedoXYZ); // pass metalSpecAlbedoRGB instead
    vec3 reflRGB = METAL_SPECULAR_REFLECTANCE(C, X, basis.nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

float evaluateMetal( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(woL.z, ior, k);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float specAlbedo = METAL_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    return f;
}

float pdfMetal( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float pdf = microfacetPDF(h, roughness) * dwh_dwo;
    return pdf;
}

float sampleMetal( in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 XYZ,
                   inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(woL.z, ior, k);
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
    if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
    float D = microfacetEval(m, roughness);
    float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
    float specAlbedo = METAL_SPEC_REFL_EVAL(X, woL, basis, XYZ);
    float f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    float dwh_dwo; // Jacobian of the half-direction mapping
    dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
    pdfOut = microfacetPDF(m, roughness) * dwh_dwo;
    return f;
}


// ****************************        Surface        ****************************

float SURFACE_DIFFUSE_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 woW, in vec3 XYZ)
{
    vec3 C = xyzToRgb(surfaceDiffuseAlbedoXYZ); // pass surfaceDiffuseAlbedoRGB instead
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(C, X, nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

float SURFACE_SPEC_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 woW, in vec3 XYZ)
{
    vec3 C = xyzToRgb(surfaceSpecAlbedoXYZ); // pass surfaceSpecAlbedoRGB instead
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(C, X, nW, woW);
    vec3 reflSPC = rgbToSpectrum(reflRGB);
    return max(dot(XYZ, reflSPC), 0.0);
}

// Fast path for non-transmissive surface
float fresnelDielectricReflectanceFast(in float cosi, in float ior)
{
    float ei = 1.0;
    float et = max(1.0, ior);
    float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));
    float cost = sqrt(max(0.0, 1.0 - sint*sint));
    float cosip = abs(cosi);
    float rParallel      = (et*cosip - ei*cost) / (et*cosip + ei*cost);
    float rPerpendicular = (ei*cosip - et*cost) / (ei*cosip + et*cost);
    return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}

float evaluateSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    float diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, XYZ);
    float    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, XYZ);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    f += diffuseAlbedo/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    vec3 woW = localToWorld(woL, basis);
    float diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, XYZ);
    float    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, XYZ);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float diffusePdf = pdfHemisphere(wiL);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    float sum = max(Fr*specAlbedo + diffuseAlbedo, DENOM_TOLERANCE);
    float specProb = Fr*specAlbedo/sum;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

float sampleSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 XYZ,
                    inout vec3 wiL, inout float pdfOut, inout vec4 rnd)
{
    vec3 woW = localToWorld(woL, basis);
    float diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, XYZ);
    float    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, XYZ);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float sum = max(specAlbedo + diffuseAlbedo, DENOM_TOLERANCE);
    float specProb = specAlbedo/sum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        wiL = sampleHemisphere(rnd, pdfOut);
        pdfOut *= (1.0-specProb);
        return diffuseAlbedo/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
        if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
        float f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = specProb * microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
}

// ****************************        BSDF common interface        ****************************

float evaluateBsdf( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ,
                    inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    evaluateSurface(X, basis, woL, wiL,                XYZ); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      evaluateMetal(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return evaluateDielectric(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
}

float sampleBsdf( in vec3 X, in Basis basis, in vec3 woL, in int material, in float wavelength_nm, in vec3 XYZ,
                  inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    sampleSurface(X, basis, woL,                XYZ, wiL, pdfOut, rnd); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      sampleMetal(X, basis, woL, wavelength_nm, XYZ, wiL, pdfOut, rnd); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return sampleDielectric(X, basis, woL, wavelength_nm, XYZ, wiL, pdfOut, rnd); }
#endif
}

float pdfBsdf( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    pdfSurface(X, basis, woL, wiL,                XYZ); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      pdfMetal(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return pdfDielectric(X, basis, woL, wiL, wavelength_nm, XYZ); }
#endif
}

#endif // HAS_GEOMETRY

#ifdef HAS_VOLUME

float VOLUME_ALBEDO_EVAL(in vec3 X, in vec3 XYZ)
{
    vec3 sigma_s = VOLUME_SCATTERING_COLOR(volumeScatteringColorRGB, X);
    vec3 sigma_a = VOLUME_ABSORPTION_COLOR(volumeAbsorptionColorRGB, X);
    vec3 albedo = sigma_s / max(sigma_s + sigma_a, 1.0e-8);
    vec3 albedoSPC = rgbToSpectrum(albedo);
    return max(dot(albedoSPC, XYZ), 0.0);
}

float VOLUME_EXTINCTION_EVAL(in vec3 X, in vec3 XYZ)
{
    vec3 sigma_s = VOLUME_SCATTERING_COLOR(volumeScatteringColorRGB, X);
    vec3 sigma_a = VOLUME_ABSORPTION_COLOR(volumeAbsorptionColorRGB, X);
    float extinction = volumeExtinction * VOLUME_EXTINCTION(volumeExtinction, X);
    vec3 sigma_t = extinction * (sigma_s + sigma_a) / lengthScale;
    vec3 sigmaSPC = rgbToSpectrum(sigma_t);
    return max(0.0, dot(sigmaSPC, XYZ));
}

float VOLUME_EXTINCTION_MAX_EVAL()
{
    return VOLUME_EXTINCTION_MAX(volumeExtinction);
}

float VOLUME_EMISSION_EVAL(in vec3 X, in vec3 XYZ)
{
    vec3 emission = VOLUME_EMISSION(volumeEmissionColorRGB * volumeEmission, X);
    vec3 emissionSPC = rgbToSpectrum(emission);
    return max(0.0, dot(emissionSPC, XYZ));
}

float VOLUME_ANISOTROPY_EVAL(in vec3 X)
{
    return VOLUME_ANISOTROPY(volumeAnisotropy, X);
}

#endif


///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit over specified segment
bool traceDistance(in vec3 start, in vec3 dir, float maxDist,
                   inout vec3 hit, inout int material)
{
    float minMarch = minLengthScale;
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
    float sdf_volum = abs(SDF_VOLUME(start));     sdf = min(sdf, sdf_volum);
#endif

    float InitialSign = sign(sdf);
    float t = InitialSign * sdf; // (always take the first step along the ray direction)
    if (t>=maxDist) return false;
    for (int n=0; n<__MAX_MARCH_STEPS__; n++)
    {
        vec3 pW = start + t*dir;
        sdf = HUGE_VAL;
#ifdef HAS_SURFACE
        sdf_surfa = abs(SDF_SURFACE(pW));    if (sdf_surfa<minMarch) { material = MAT_SURFA; hit = start + t*dir; return true; } sdf = min(sdf, sdf_surfa);
#endif
#ifdef HAS_METAL
        sdf_metal = abs(SDF_METAL(pW));      if (sdf_metal<minMarch) { material = MAT_METAL; hit = start + t*dir; return true; } sdf = min(sdf, sdf_metal);
#endif
#ifdef HAS_DIELECTRIC
        sdf_diele = abs(SDF_DIELECTRIC(pW)); if (sdf_diele<minMarch) { material = MAT_DIELE; hit = start + t*dir; return true; } sdf = min(sdf, sdf_diele);
#endif
#ifdef HAS_VOLUME
        sdf_volum = abs(SDF_VOLUME(pW));     if (sdf_volum<minMarch) { material = MAT_VOLUM; hit = start + t*dir; return true; } sdf = min(sdf, sdf_volum);
#endif
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf;
        if (t>=maxDist) return false;
    }
    return !maxStepsIsMiss;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir,
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_INVAL;
    return traceDistance(start, dir, maxMarchDist, hit, material);
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

// MC-estimates the amount of light transmitted along an infinite ray
float Visibility(in vec3 pW, in vec3 rayDir, in vec3 XYZ, inout vec4 rnd, bool inVolume)
{
    vec3 pW_next;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
#ifndef HAS_VOLUME
    return float(!hit);
#else
    // Deal with all volume interactions until the ray either hits a surface or escapes
    const int maxVolHits = 16;
    int volHits = 0;
    float sigma_t_max = VOLUME_EXTINCTION_MAX_EVAL();
    float inv_sigma_t_max = 1.0/sigma_t_max;
    float eps = 3.0*minLengthScale;
    while ((inVolume || hitMaterial==MAT_VOLUM) && volHits<maxVolHits)
    {
        // If ray lies in the volume interior
        if (inVolume)
        {
            // Do Woodcock tracking to check for a possible volume scattering event along the in-volume segment pW -> pW_next
            // (except if the current ray lies inside a dielectric, in which case there is no volumetric scattering, i.e. the dielectric displaces the volume)
            vec3 pScatter;
            float sigma_t;
            bool interacted = false;
            float segmentLength = length(pW_next - pW);
            float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max;
            int steps = 0;
            while (distanceMarched<segmentLength-eps && steps<__MAX_VOLUME_STEPS__)
            {
                pScatter = pW + distanceMarched*rayDir;
                sigma_t = VOLUME_EXTINCTION_EVAL(pScatter, XYZ);
                interacted = bool(rand(rnd) < sigma_t * inv_sigma_t_max);
                if (interacted) break;
                distanceMarched += -log(rand(rnd)) * inv_sigma_t_max;
                steps++;
            }

            if (interacted) return 0.0; // on absorption/scattering, return zero visibility
            else
            {
                pW = pW_next;
                if (hitMaterial==MAT_VOLUM) // if hitMaterial is volume boundary, displace pW into volume exterior
                {
                    vec3 nW = normal(pW, MAT_VOLUM);
                    pW += nW * eps;
                    inVolume = false;
                }
                else // if hitMaterial is a surface, return zero visibility
                    return 0.0;
            }
        }

        // else if ray passed through the volume exterior and hit the volume boundary
        else // !inVolume
        {
            pW = pW_next; // advance to the volume boundary
            vec3 nW = normal(pW, MAT_VOLUM); // displace pW into volume interior
            pW -= nW * eps;
            inVolume = true;
        }

        // Re-trace from new pW
        hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
        volHits++;
    }
    return float(!hit || volHits==maxVolHits); // if no hit or scattering, return visibility 1.0 (transparent)
#endif // HAS_VOLUME
}


////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dir)
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


float environmentRadiance(in vec3 dir, in vec3 XYZ)
{
    vec3 RGB_sky = environmentRadianceRGB(dir);
    vec3 XYZ_spec = rgbToSpectrum(RGB_sky); // convert to radiance at the given wavelength
    return dot(XYZ, XYZ_spec);
}

vec3 sampleSunDir(inout vec4 rnd)
{
    float costhetamax = cos(sunAngularSize * M_PI/180.0);
    float costheta = 1.0 - rand(rnd)*(1.0-costhetamax);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    return localToWorld(vec3(x, y, z), sunBasis);
}

float sunRadiance(in vec3 dir, in vec3 XYZ)
{
    if (dot(dir, sunDir) < cos(sunAngularSize*M_PI/180.0)) return 0.0;
    vec3 RGB_sun = sunPower * sunColor;
    vec3 XYZ_sun = rgbToXyz(RGB_sun); // @todo: just supply XYZ_sun
    vec3 XYZ_spec = xyz_to_spectrum(XYZ_sun); // convert to radiance at the given wavelength
    return dot(XYZ, XYZ_spec);
}

#ifdef HAS_GEOMETRY

float sampleLightAtSurface(Basis basis, in vec3 XYZ, inout vec4 rnd, inout vec3 wiL, inout vec3 wiW, inout float lightPdf)
{
    // Light sampling (choose either sun or sky)
    vec3 wiW_sun = sampleSunDir(rnd);
    vec3 wiL_sun = worldToLocal(wiW_sun, basis);
    float sunPdf = 1.0; // convenient to decouple total sun power from its angular size
    float sunWeight = sunPower * max(0.0, wiL_sun.z);

    float skyPdf;
    vec3 wiL_sky = sampleHemisphere(rnd, skyPdf);
    float skyWeight = skyPower * max(0.0, wiL_sky.z);

    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    bool chooseSun = (rand(rnd) <= sunProb);
    if (chooseSun)
    {
        lightPdf = sunPdf * sunProb;
        wiL = wiL_sun;
        wiW = wiW_sun;
        return sunRadiance(wiW_sun, XYZ);
    }
    lightPdf = skyPdf * max(PDF_EPSILON, 1.0-sunProb);
    wiL = wiL_sky;
    wiW = localToWorld(wiL_sky, basis);
    return environmentRadiance(wiW, XYZ);
}

// Estimate direct radiance at the given surface vertex
float directSurfaceLighting(in vec3 pW, Basis basis, in vec3 woW, in int material,
                            float wavelength_nm, in vec3 XYZ, inout vec4 rnd, bool inVolume, inout float lightPdf)
{
    vec3 wiL, wiW; // direction of sampled direct light (*towards* the light)
    float Li = sampleLightAtSurface(basis, XYZ, rnd, wiL, wiW, lightPdf);
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    float V = Visibility(pW+dPw, wiW, XYZ, rnd, inVolume);
    Li *= abs(1.0 - shadowStrength*(1.0-V));

    // Apply MIS weight with the BSDF pdf for the sampled direction
    vec3 woL = worldToLocal(woW, basis);
    float bsdfPdf = pdfBsdf(pW, basis, woL, wiL, material, wavelength_nm, XYZ);
    if (bsdfPdf<PDF_EPSILON) return 0.0;
    float f = evaluateBsdf(pW, basis, woL, wiL, material, wavelength_nm, XYZ, rnd);
    float misWeight = powerHeuristic(lightPdf, bsdfPdf);
    float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li * abs(dot(wiW, basis.nW)) * misWeight;
}

#endif // HAS_GEOMETRY

#ifdef HAS_VOLUME

float sampleLightInVolume(in vec3 XYZ, inout vec4 rnd, inout vec3 wiW, inout float lightPdf)
{
    // Light sampling (choose either sun or sky)
    vec3 wiW_sun = sampleSunDir(rnd);
    float sunPdf = 1.0; // convenient to decouple total sun power from its angular size
    float sunWeight = sunPower * sunPdf;
    float skyPdf;
    vec3 wiW_sky = sampleSphere(rnd, skyPdf);
    float skyWeight = skyPower * skyPdf;
    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    bool chooseSun = (rand(rnd) <= sunProb);
    if (chooseSun)
    {
        lightPdf = sunPdf * sunProb;
        wiW = wiW_sun;
        return sunRadiance(wiW_sun, XYZ);
    }
    lightPdf = skyPdf * max(PDF_EPSILON, 1.0-sunProb);
    wiW = wiW_sky;
    return environmentRadiance(wiW, XYZ);
}

// Estimate direct radiance at the given volumetric vertex
float directVolumeLighting(in vec3 pW, in vec3 woW, in vec3 XYZ, inout vec4 rnd, bool inVolume)
{
    vec3 wiW; // direction of sampled direct light (*towards* the light)
    float lightPdf;
    float Li = sampleLightInVolume(XYZ, rnd, wiW, lightPdf);
    float V = Visibility(pW, wiW, XYZ, rnd, inVolume);
    Li *= abs(1.0 - shadowStrength*(1.0-V));
    float f = phaseFunction(dot(woW, -wiW));
    float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li;
}

#endif // HAS_VOLUME

////////////////////////////////////////////////////////////////////////////////
// Pathtracing integrator
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
    float theta = 2.0*M_PI*rand(rnd);
    vec3 lensPos = camPos + lensRadial*(-camX*cos(theta) + camY*sin(theta));
    primaryStart = lensPos;
    primaryDir = normalize(focalPlaneHit - lensPos);
}

void pathtrace(vec2 pixel, vec4 rnd) // the current pixel
{
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

    // Sample photon wavelength via the inverse CDF of the emission spectrum.
    // We limit the sampled wavelengths to between about 420nm and 700nm to avoid color mismatch at low sample count
    float xi = 0.1+0.8*(0.5+floor(float(wavelengthSamples)*rand(rnd))) / float(wavelengthSamples);
    float w = texture(ICDF, vec2(xi, 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;

    // Convert wavelength to XYZ tristimulus
    vec3 XYZ = texture(WavelengthToXYZ, vec2(w, 0.5)).rgb;
    vec3 RGB = clamp(xyzToRgb(XYZ), 0.0, 1.0);
    
    // Jitter over pixel
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    // Compute world ray direction for this fragment
    vec3 primaryStart, primaryDir;
    constructPrimaryRay(pixel, rnd, primaryStart, primaryDir);

    // Perform pathtrace to estimate the primary ray radiance, L
    float L = 0.0;
    float throughput = 1.0;
    float misWeight = 1.0; // For MIS book-keeping
    vec3 pW = primaryStart;
    vec3 rayDir = primaryDir; // (opposite to light direction)
    bool inVolume = false;
#ifdef HAS_VOLUME
    float sigma_t_max = VOLUME_EXTINCTION_MAX_EVAL();
    float inv_sigma_t_max = 1.0/sigma_t_max;
    inVolume = SDF_VOLUME(primaryStart) < 0.0;
#endif
    float eps = 3.0*minLengthScale;
#ifdef HAS_DIELECTRIC
    bool inDielectric = SDF_DIELECTRIC(primaryStart) < 0.0;
#endif

    for (int vertex=0; vertex<=__MAX_BOUNCES__; ++vertex)
    {
        // Raycast along current propagation direction rayDir, from current vertex pW to pW_next
        vec3 pW_next;
        int hitMaterial;
        bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
        float rayLength = maxLengthScale;
        if (hit)
        {
            rayLength = length(pW_next - pW);
        }

#ifdef HAS_VOLUME
        // Deal with all volume interactions until the ray either hits a surface or escapes
        int volHits = 0;
        bool absorb = false;
        while ((inVolume || hitMaterial==MAT_VOLUM) && volHits<(1+__MAX_SCATTERS__))
        {
            // If ray lies in the volume interior
            if (inVolume)
            {
                float segmentLength = length(pW_next - pW);

                // Do Woodcock tracking to check for a possible volume absorption/scattering event along the in-volume segment pW -> pW_next
                // (except if the current ray lies inside a dielectric, in which case there is no volumetric scattering, i.e. the dielectric displaces the volume)
                vec3 pScatter;
                float sigma_t;
                bool interacted = false;
#ifdef HAS_DIELECTRIC
                if (!inDielectric)
#endif
                {
                    float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max;
                    int steps = 0;
                    while (distanceMarched<segmentLength-eps && steps<__MAX_VOLUME_STEPS__)
                    {
                        pScatter = pW + distanceMarched*rayDir;
                        sigma_t = VOLUME_EXTINCTION_EVAL(pScatter, XYZ);
                        interacted = bool(rand(rnd) < sigma_t * inv_sigma_t_max);
                        if (interacted) break;
                        distanceMarched += -log(rand(rnd)) * inv_sigma_t_max;
                        steps++;
                    }
                }

                // Absorption/scattering event
                if (interacted)
                {
                    vec3 wiW = samplePhaseFunction(rayDir, rnd); // sample scattered dir (NB, PF is the angle PDF so the MC PDF denom. cancels it in the L estimator)
                    float emission = VOLUME_EMISSION_EVAL(pScatter, XYZ);
                    float albedo = VOLUME_ALBEDO_EVAL(pScatter, XYZ);
                    absorb = (rand(rnd) > albedo);
                    if (absorb)
                    {
                        L += throughput * emission; // add volume emission term
                        break;
                    }

                    L += throughput * directVolumeLighting(pScatter, -rayDir, XYZ, rnd, inVolume); // add contribution due to direct lighting at vertex
                    throughput *= albedo; // update throughput due to scattering
                    if (throughput < THROUGHPUT_EPSILON) break;
                    pW = pScatter; // continue to next vertex
                    rayDir = wiW;
                }

                // No scattering
                else
                {
                    pW = pW_next;
                    if (hitMaterial==MAT_VOLUM) // if hitMaterial is volume boundary, displace pW into volume exterior
                    {
                        vec3 nW = normal(pW, MAT_VOLUM);
                        pW += nW * eps;
                        inVolume = false;
                    }
                    else        // if hitMaterial is a surface, break and handle it (note, still inside volume)
                        break;  // (or if !hit, that's a numerical glitch as we should be inside a closed volume boundary)
                }
            }

            // else if ray passed through the volume exterior and hit the volume boundary
            else // !inVolume
            {
                pW = pW_next; // advance to the volume boundary
                vec3 nW = normal(pW, MAT_VOLUM); // displace pW into volume interior
                pW -= nW * eps;
                inVolume = true;
            }

            // Re-trace from new pW
            hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
            volHits++;
        }

        // If ray was absorbed by the volume, terminate
        if (absorb) break;
        if (volHits == __MAX_SCATTERS__) hit = false;
#endif // HAS_VOLUME

        // If this ray didn't scatter and missed all geometry, add environment light term and terminate path
        if (!hit)
        {
            float Li = 0.0;
            bool sunVisible = true; // @todo: parameter
            if (!(vertex==0 && !envMapVisible))      Li += environmentRadiance(rayDir, XYZ);
            if (!(vertex==0 && !sunVisibleDirectly)) Li += sunRadiance(rayDir, XYZ);
            L += throughput * Li * misWeight;
            break;
        }

#ifndef HAS_GEOMETRY
        break;
#else
        // If the current ray lies inside a dielectric, apply Beer's law for absorption
#ifdef HAS_DIELECTRIC
        if (inDielectric)
        {
            throughput *= exp(-rayLength*dot(dieleAbsorptionRGB, RGB));
        }
#endif

        // This ray didn't scatter but hit some geometry, so deal with the surface interaction.
        // First, compute normal at current surface vertex
        pW = pW_next;
        vec3 nW = normal(pW, hitMaterial);
        Basis basis = makeBasis(nW);

        // Sample BSDF for the next bounce direction
        vec3 woW = -rayDir;
        vec3 woL = worldToLocal(woW, basis);
        vec3 wiL;
        float bsdfPdf;
        float f = sampleBsdf(pW, basis, woL, hitMaterial, wavelength_nm, XYZ, wiL, bsdfPdf, rnd);
        vec3 wiW = localToWorld(wiL, basis);
        rayDir = wiW; // Update ray direction

        // Detect dielectric transmission
#ifdef HAS_DIELECTRIC
        if (hitMaterial==MAT_DIELE && dot(woW, nW)*dot(wiW, nW) < 0.0)
        {
            inDielectric = !inDielectric;
        }
#endif

        // Update path throughput
        float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, bsdfPdf));
        throughput *= fOverPdf * abs(dot(wiW, nW));
        if (throughput < THROUGHPUT_EPSILON) break;

        // Add direct lighting term at current surface vertex
        float lightPdf = 0.0;
#ifdef HAS_DIELECTRIC
        if (!inDielectric)
#endif
        {
            L += throughput * directSurfaceLighting(pW, basis, woW, hitMaterial, wavelength_nm, XYZ, rnd, inVolume, lightPdf);
        }

        // Prepare for tracing the bounce ray
        pW += nW * sign(dot(wiW, nW)) * eps;           // perturb vertex into half-space of scattered ray
        misWeight = powerHeuristic(bsdfPdf, lightPdf); // compute MIS weight for bounce ray
#endif // HAS_GEOMETRY
    }

    // Compute tristimulus contribution from estimated radiance
    vec3 colorXYZ = XYZ * L;

    // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;

    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}

void main()
{
    INIT();
    vec4 rnd = texture(RngData, vTexCoord);
    pathtrace(gl_FragCoord.xy, rnd);
}
`,

'pathtracer-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'pick-fragment-shader': `#version 300 es
precision highp float;

layout(location = 0) out vec4 gbuf_pick;

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
uniform bool maxStepsIsMiss;
uniform vec2 mousePick;

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

////////////////////////////////////////////////////////////////////////////////
// Picking program: computes first hit distance based on mouse location
////////////////////////////////////////////////////////////////////////////////

void constructPickRay(in vec2 pixel,
                      inout vec3 primaryStart, inout vec3 primaryDir)
{
    // Compute world ray direction for given (possibly jittered) fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
    float fh = tan(0.5*radians(camFovy)); // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    primaryDir = normalize(camDir + s);
    primaryStart = camPos;
}

void main()
{
    INIT();

    vec2 pixel = mousePick;
    vec3 primaryStart, primaryDir;
    constructPickRay(pixel, primaryStart, primaryDir);

     // Raycast to first hit point
    vec3 pW;
    int hitMaterial;
    bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxLengthScale);
    float d = maxLengthScale;
    if (hit)
    {
        d = length(pW - primaryStart);
    }

    gbuf_pick = vec4(d, hitMaterial, 0, 0);
}
`,

'pick-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'tonemapper-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Radiance;
in vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float whitepoint;

out vec4 outputColor;


void constrain_rgb(inout vec3 RGB)
{
    float w;
    w = (0.0 < RGB.r) ? 0.0 : RGB.r;
    w = (w   < RGB.g) ? 0.0 : RGB.g;
    w = (w   < RGB.b) ? 0.0 : RGB.b;
    w = -w;
    if (w>0.0)
    {
        RGB.r += w; RGB.g += w; RGB.b += w;
    }
}

void main()
{
    vec3 L = exposure * texture(Radiance, vTexCoord).rgb;
    float X = L.x;
    float Y = L.y;
    float Z = L.z;
    float sum = X + Y + Z;
    float x = X / sum;
    float y = Y / sum; 

    // compute Reinhard tonemapping scale factor
    float scale = (1.0 + Y/(whitepoint*whitepoint)) / (1.0 + Y);
    Y *= scale;
    X = x * Y / y;
    Z = (1.0 - x - y) * (Y / y);

    // convert XYZ tristimulus to sRGB color space
    vec3 RGB;
    RGB.r =  3.2406*X - 1.5372*Y - 0.4986*Z;
    RGB.g = -0.9689*X + 1.8758*Y + 0.0415*Z;
    RGB.b =  0.0557*X - 0.2040*Y + 1.0570*Z;

    // deal with out-of-gamut RGB.
    constrain_rgb(RGB);

    // apply gamma correction
    vec3 S = pow(abs(RGB), vec3(invGamma));

    outputColor = vec4(S, 0.0);
}
`,

'tonemapper-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;
out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

}