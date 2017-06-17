
#extension GL_EXT_draw_buffers : require
precision highp float;

uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D RadianceBlocks;   // 4
uniform sampler2D iorTex;           // 5 (for metals)
uniform sampler2D kTex;             // 6 (for metals)
uniform sampler2D envMap;           // 7
varying vec2 vTexCoord;

uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camNear;
uniform float camFar;
uniform float camFovy; // degrees 
uniform float camZoom;
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

uniform float minScale;
uniform float maxScale;
uniform float skyPower;
uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapRotation;
uniform float gamma;
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

#define MAT_VACUU  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2

#define M_PI 3.1415926535897932384626433832795

//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

SHADER

///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit over specified segment
bool traceDistance(in vec3 start, in vec3 dir, float maxDist,
                   inout vec3 hit, inout int material)
{
    float minMarchDist = minScale;
    float sdf_diele = abs(SDF_DIELECTRIC(start));
    float sdf_metal = abs(SDF_METAL(start));
    float sdf_surfa = abs(SDF_SURFACE(start));
    float sdf = min(sdf_diele, min(sdf_metal, sdf_surfa));
    float InitialSign = sign(sdf);

    float t = 0.0;  
    int iters=0;
    for (int n=0; n<MAX_MARCH_STEPS; n++)
    {
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root. 
        t += InitialSign * sdf;
        if (t>=maxDist) break;
        vec3 pW = start + t*dir;    
        sdf_diele = abs(SDF_DIELECTRIC(pW)); if (sdf_diele<minMarchDist) { material = MAT_DIELE; break; }
        sdf_metal = abs(SDF_METAL(pW));      if (sdf_metal<minMarchDist) { material = MAT_METAL; break; }
        sdf_surfa = abs(SDF_SURFACE(pW));    if (sdf_surfa<minMarchDist) { material = MAT_SURFA; break; }
        sdf = min(sdf_diele, min(sdf_metal, sdf_surfa));
        iters++;
    }
    hit = start + t*dir;
    if (t>=maxDist) return false;
    if (maxStepsIsMiss && iters>=MAX_MARCH_STEPS) return false;
    return true;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir, 
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_VACUU;
    return traceDistance(start, dir, maxMarchDist, hit, material);
}

// (whether occluded along infinite ray)
bool Occluded(in vec3 start, in vec3 dir)
{
    float eps = 3.0*minScale;
    vec3 delta = eps * dir;
    int material;
    vec3 hit;
    return traceRay(start+delta, dir, hit, material, maxScale);
}

vec3 normal(in vec3 pW, int material)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = 2.0*minScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 xyyp = pW+e.xyy; vec3 xyyn = pW-e.xyy;
    vec3 yxyp = pW+e.yxy; vec3 yxyn = pW-e.yxy;
    vec3 yyxp = pW+e.yyx; vec3 yyxn = pW-e.yyx;
    vec3 N;
    if      (material==MAT_DIELE) { N = vec3(SDF_DIELECTRIC(xyyp)-SDF_DIELECTRIC(xyyn), SDF_DIELECTRIC(yxyp)-SDF_DIELECTRIC(yxyn), SDF_DIELECTRIC(yyxp) - SDF_DIELECTRIC(yyxn)); }
    else if (material==MAT_METAL) { N = vec3(SDF_METAL(xyyp)     -SDF_METAL(xyyn),      SDF_METAL(yxyp)     -SDF_METAL(yxyn),      SDF_METAL(yyxp)      - SDF_METAL(yyxn)); }
    else                          { N = vec3(SDF_SURFACE(xyyp)   -SDF_SURFACE(xyyn),    SDF_SURFACE(yxyp)   -SDF_SURFACE(yxyn),    SDF_SURFACE(yyxp)    - SDF_SURFACE(yyxn)); }
    return normalize(N);
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
        vec3 RGB = texture2D(envMap, vec2(u,v)).rgb;
        RGB.r = pow(RGB.r, gamma);
        RGB.g = pow(RGB.g, gamma);
        RGB.b = pow(RGB.b, gamma);
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

void main()
{
    INIT();

    vec4 rnd = texture2D(RngData, vTexCoord);
    vec2 pixel = gl_FragCoord.xy;
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    vec3 primaryStart, primaryDir;
    constructPrimaryRay(pixel, rnd, primaryStart, primaryDir);

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int hitMaterial;
    bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxScale);
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
    vec4 oldL = texture2D(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;

    gl_FragData[0] = vec4(newL, newN);
    gl_FragData[1] = rnd;
}

