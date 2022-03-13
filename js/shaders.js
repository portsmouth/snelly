var Shaders = {

'ao-fragment-shader-backup': `#version 300 es
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

// Camera parameters
uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

// Rendering  parameters
uniform float radianceClamp;
uniform float skipProbability;
uniform bool maxStepsIsMiss;

// Length scales
uniform float lengthScale;                      // UPLOAD!
uniform float minLengthScale;
uniform float maxLengthScale;

// Surface material parameters
uniform vec3 surfaceDiffuseAlbedoRGB;

// Sky parameters
uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapPhiRotation;
uniform float envMapThetaRotation;
uniform float envMapTransitionAngle;
uniform float skyPower;
uniform vec3 skyTintUp;
uniform vec3 skyTintDown;

// Sun parameters
uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;                      // UPLOAD!

// Atmosphere constants
uniform vec3 atmosphereBoundsMin;                      // UPLOAD!
uniform vec3 atmosphereBoundsMax;                      // UPLOAD!
uniform vec3 atmosphereExtinctionCoeffRGB;                      // UPLOAD!
uniform vec3 atmosphereScatteringCoeffRGB;                      // UPLOAD!
uniform float atmosphereMFP;                      // UPLOAD!
uniform float atmosphereAnisotropy;                      // UPLOAD!

#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5
#define RADIANCE_EPSILON 1.0e-5

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2

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
float Transmittance(in vec3 start, in vec3 dir)
{
    int material;
    vec3 hit;
    return traceRay(start, dir, hit, material, maxLengthScale) ? 0.0 : 1.0;
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

float averageComponent(in vec3 v)
{
    return (v.r + v.g + v.b)/3.0;
}


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphereCosineWeighted(inout vec4 rnd, inout float pdf)
{
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}

#ifdef HAS_ATMOSPHERE
vec3 sampleSphereUniformly(inout vec4 rnd)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    // pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}
#endif

vec3 SURFACE_DIFFUSE_REFL_RGB(in vec3 X, in vec3 nW, in vec3 woW)
{
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, woW);
    return reflRGB;
}

#ifdef HAS_SURFACE_NORMALMAP
vec3 perturbNormal(in vec3 X, in Basis basis, int material)
{
    if (material==MAT_SURFA)
    {
        vec3 nL = SURFACE_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL), basis);
    }
    return basis.nW;
}
#endif

//////////////////////////////////////////////
// Sky
//////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dirW)
{
    float rot_phi   = M_PI*envMapPhiRotation/180.0;
    float rot_theta = M_PI*envMapThetaRotation/180.0;
    vec3 sky_pole = vec3(sin(rot_theta), cos(rot_theta), 0.0);
    Basis sky_basis = makeBasis(sky_pole);
    vec3 dirL = worldToLocal(dirW, sky_basis);
    float phi = atan(dirL.y, dirL.x) + M_PI + rot_phi;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = acos(dirL.z);
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 RGB = vec3(1.0);
    if (haveEnvMap)
        RGB = texture(envMap, vec2(u,v)).rgb;
    float t = dot(dirW, sky_pole);
    float tt = envMapTransitionAngle/180.0;
    RGB *= skyPower * mix(skyTintDown, skyTintUp, smoothstep(-tt, tt, t));
    return RGB;
}

vec3 environmentRadiance(in vec3 dir)
{
    vec3 RGB_sky = environmentRadianceRGB(dir);
    return RGB_sky;
}

vec3 sampleSkyAtSurface(in Basis basis, inout vec4 rnd,
                        inout vec3 woutputW, inout float pdfDir)
{
    vec3 woutputL = sampleHemisphereCosineWeighted(rnd, pdfDir);
    woutputW = localToWorld(woutputL, basis);
    return environmentRadiance(woutputW);
}

#ifdef HAS_ATMOSPHERE
vec3 sampleSkyInVolume(inout vec4 rnd,
                       inout vec3 woutputW, inout float pdfDir)
{
    if (skyPower<RADIANCE_EPSILON)
        return vec3(0.0);
    woutputW = sampleSphereUniformly(rnd);
    pdfDir = 1.0/(4.0*M_PI);
    return environmentRadiance(woutputW);
}
#endif

//////////////////////////////////////////////
// Sun
//////////////////////////////////////////////

Basis sunBasis;

vec3 sampleSunDir(inout vec4 rnd, inout float pdfDir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    float theta = theta_max * sqrt(rand(rnd));
    float costheta = cos(theta);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    pdfDir = 1.0/solid_angle;
    return localToWorld(vec3(x, y, z), sunBasis);
}

float pdfSun(in vec3 dir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return 0.0;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    return 1.0/solid_angle;
}

vec3 sunRadiance(in vec3 dir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return vec3(0.0);
    vec3 RGB_sun = sunPower * sunColor;
    return RGB_sun;
}

vec3 sampleSunAtSurface(in Basis basis, inout vec4 rnd,
                        inout vec3 woutputW, inout float pdfDir)
{
    woutputW = sampleSunDir(rnd, pdfDir);
    vec3 woutputL = worldToLocal(woutputW, basis);
    if (woutputL.z < 0.0) return vec3(0.0);
    return sunRadiance(woutputW);
}

#ifdef HAS_ATMOSPHERE
vec3 sampleSunInVolume(inout vec4 rnd,
                       inout vec3 woutputW, inout float pdfDir)
{
    if (sunPower<RADIANCE_EPSILON)
        return vec3(0.0);
    woutputW = sampleSunDir(rnd, pdfDir);
    return sunRadiance(woutputW);
}
#endif


//////////////////////////////////////////////
// Atmosphere
//////////////////////////////////////////////

#ifdef HAS_ATMOSPHERE

vec3 VOLUME_SCATTERING_EVAL()
{
    return atmosphereScatteringCoeffRGB/(atmosphereMFP*lengthScale);
}

vec3 VOLUME_EXTINCTION_EVAL()
{
    return atmosphereExtinctionCoeffRGB/(atmosphereMFP*lengthScale);
}

#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }

bool atmosphereHit( in vec3 rayPos, in vec3 rayDir,
                    inout float t0, inout float t1 )
{
    vec3 dL = 1.0/rayDir;
    vec3 lo = (atmosphereBoundsMin*lengthScale - rayPos) * dL;
    vec3 hi = (atmosphereBoundsMax*lengthScale - rayPos) * dL;
    sort2(lo, hi);
    bool hit = !( lo.x>hi.y || lo.y>hi.x || lo.x>hi.z || lo.z>hi.x || lo.y>hi.z || lo.z>hi.y );
    t0 = max(max(lo.x, lo.y), lo.z);
    t1 = min(min(hi.x, hi.y), hi.z);
    return hit;
}

bool atmosphereSegment(in vec3 pW, in vec3 rayDir, float segmentLength, // input segment
                       inout float t0, inout float t1)                  // portion of input segment within atmosphere
{
    if ( !atmosphereHit(pW, rayDir, t0, t1) )
        return false;
    t0 = max(0.0, t0);
    if (t1<0.0 || segmentLength<t0)
        return false;
    t1 = min(segmentLength, t1);
    return true;
}

// Return the amount of light transmitted (per-channel) along a known "free" segment (i.e. free of geometry)
vec3 transmittanceOverFreeSegment(in vec3 pW, in vec3 rayDir, float segmentLength)
{
    vec3 extinction = VOLUME_EXTINCTION_EVAL();
    if (length(extinction) == 0.f)
        return vec3(1.0);
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return vec3(1.0);
    vec3 opticalDepth = (t1 - t0) * extinction;
    vec3 Tr = exp(-opticalDepth);
    return Tr;
}

// Return the amount of light transmitted (per-channel) along a given segment
// (zero if occluded by geometry).
vec3 transmittanceOverSegment(in vec3 pW, in vec3 rayDir, float segmentLength)
{
    vec3 pW_surface;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_surface, hitMaterial, segmentLength);
    if (hit)
        return vec3(0.0);
    else
        return transmittanceOverFreeSegment(pW, rayDir, segmentLength);
}

float phaseFunction(float mu, float anisotropy)
{
    float g = anisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}

vec3 directVolumeLighting(in vec3 pW, in vec3 rayDir, in vec3 rgb, inout vec4 rnd)
{
    vec3 Ldirect = vec3(0.0);
    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        vec3 woutputW;
        float skyPdf;
        vec3 Li = sampleSkyInVolume(rnd, woutputW, skyPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            vec3 Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale);
            Li *= Tr;
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, skyPdf);
            }
        }
    }
    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        vec3 woutputW;
        float sunPdf;
        vec3 Li = sampleSunInVolume(rgb, rnd, woutputW, sunPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            vec3 Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale);
            Li *= Tr;
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, sunPdf);
            }
        }
    }
    return min(vec3(radianceClamp), Ldirect);
}

vec3 atmosphericInscatteringRadiance(in vec3 pW, in vec3 rayDir, in float segmentLength, inout vec4 rnd)
{
    vec3 Ls = vec3(0.0);
    vec3 extinction = VOLUME_EXTINCTION_EVAL();
    vec3 scattering = VOLUME_SCATTERING_EVAL();
    float extinction_norm = averageComponent(extinction);
    float scattering_norm = averageComponent(scattering);
    if ( extinction_norm < RADIANCE_EPSILON ||
         scattering_norm < RADIANCE_EPSILON ) return Ls;

    // Find sub-segment of supplied segment over which (homogeneous) atmosphere exists
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return Ls;

    // Sample a scattering point in [t0, t1] from a PDF proportional to the transmittance (normalized over [t0, t1])
    float T01 = exp(-(t1-t0)*extinction_norm);
    float t_scatter = t0 - log(1.0 - rand(rnd)*(1.0 - T01)) / extinction_norm; // sampled scatter distance
    vec3 pW_scatter = pW + t_scatter*rayDir;                                   // sampled scatter point
    vec3 opticalDepth_scatter = (t_scatter - t0) * extinction;         // optical depth from pW -> pW_scatter
    vec3 Tr_pW = exp(-opticalDepth_scatter);                           // transmittance from pW -> pW_scatter
    float invDistancePdf = (1.0 - T01) * exp(averageComponent(opticalDepth_scatter)) / extinction_norm; // PDF of scatter distance

    // Direct lighting
    vec3 Li = directVolumeLighting(pW_scatter, rayDir, rnd);
    Ls += Tr_pW * scattering * Li * invDistancePdf; // final estimator for scattered radiance

    return Ls;
}

#endif // HAS_ATMOSPHERE


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
    vec4 rnd = texture(RngData, vTexCoord);
#ifdef INTERACTIVE_MODE
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
#endif

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    INIT();

    vec2 pixel = gl_FragCoord.xy;
    vec3 RGB = vec3(0.0);
    for (int n=0; n<__MAX_SAMPLES_PER_FRAME__; ++n)
    {
        // Jitter over pixel
        vec2 pixelj = pixel + (-0.5 + vec2(rand(rnd), rand(rnd)));
        vec3 primaryStart, primaryDir;
#ifdef HAS_CUSTOM_CAMERA
        CONSTRUCT_PRIMARY_RAY(pixelj, rnd, primaryStart, primaryDir);
#else
        constructPrimaryRay(pixelj, rnd, primaryStart, primaryDir);
#endif

        // Raycast to first hit point
        vec3 pW;
        vec3 woW = -primaryDir;
        int hitMaterial;
        bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxLengthScale);

        // Closest-hit logic
        if (hit)
        {
            // Compute normal at hit point
            vec3 nW = normal(pW, hitMaterial);
            Basis basis = makeBasis(nW);
#ifdef HAS_SURFACE_NORMALMAP
            nW = perturbNormal(pW, basis, hitMaterial);
            basis = makeBasis(nW);
#endif

#ifdef HAS_ATMOSPHERE
            // Add term for single-scattering in the homogeneous atmosphere over the segment up to the closest hit
            float rayLength = length(pW - primaryStart);
            vec3 Lscatter = atmosphericInscatteringRadiance(primaryStart, primaryDir, rayLength, rnd);
            RGB += Lscatter;
#endif

            // Compute diffuse BSDF
            vec3 f = SURFACE_DIFFUSE_REFL_RGB(pW, nW, woW) / M_PI;

            // Compute direct lighting
            vec3 Ldirect = vec3(0.0);
            vec3 dPw = 3.0*minLengthScale*nW;
            if (skyPower > RADIANCE_EPSILON)
            {
                // Sky
                float skyPdf;
                vec3 woutputW;
                vec3 Li = sampleSkyAtSurface(basis, rnd, woutputW, skyPdf);
#ifdef HAS_ATMOSPHERE
                vec3 Tr = transmittanceOverSegment(pW+dPw, woutputW, maxLengthScale, rgb);
#else
                vec3 Tr = vec3(1.0);
#endif
                Ldirect += f * Tr * Li / max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, nW));
            }
            if (sunPower > RADIANCE_EPSILON)
            {
                // Sun
                float sunPdf;
                vec3 woutputW;
                vec3 Li = sampleSunAtSurface(basis, rnd, woutputW, sunPdf);
#ifdef HAS_ATMOSPHERE
                vec3 Tr = transmittanceOverSegment(pW+dPw, woutputW, maxLengthScale, rgb);
#else
                vec3 Tr = vec3(1.0);
#endif
                if (averageComponent(Li) > RADIANCE_EPSILON)
                {
                    Li *= Transmittance(pW+dPw, woutputW);
                    Ldirect += f * Tr * Li / max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, nW));
                }
            }
            RGB += Ldirect;
        }

        // Miss logic
        else
        {
            vec3 Tr = vec3(1.0);
#ifdef HAS_ATMOSPHERE
            // Add term for single-inscattering in the homogeneous atmosphere over the segment to infinity
            vec3 Lscatter = atmosphericInscatteringRadiance(primaryStart, primaryDir, maxLengthScale, rgb, rnd);
            RGB += Lscatter;
            Tr = transmittanceOverFreeSegment(primaryStart, primaryDir, maxLengthScale, rgb);
#endif
            RGB += Tr * environmentRadiance(primaryDir);
            //if (sunVisibleDirectly)
            //    RGB += Tr * sunRadiance(primaryDir);
        }

    } // sample n of __MAX_SAMPLES_PER_FRAME__

    RGB /= float(__MAX_SAMPLES_PER_FRAME__);

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

'ao-fragment-shader': `#version 300 es
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

// Camera parameters
uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

// Rendering  parameters
uniform float radianceClamp;
uniform float skipProbability;
uniform bool maxStepsIsMiss;

// Length scales
uniform float lengthScale;                      // UPLOAD!
uniform float minLengthScale;
uniform float maxLengthScale;

// Surface material parameters
uniform vec3 surfaceDiffuseAlbedoRGB;

// Sky parameters
uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapPhiRotation;
uniform float envMapThetaRotation;
uniform float envMapTransitionAngle;
uniform float skyPower;
uniform vec3 skyTintUp;
uniform vec3 skyTintDown;

// Sun parameters
uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;

// Atmosphere constants
uniform vec3 atmosphereBoundsMin;
uniform vec3 atmosphereBoundsMax;
uniform vec3 atmosphereExtinctionCoeffRGB;
uniform vec3 atmosphereScatteringCoeffRGB;
uniform float atmosphereMFP;
uniform float atmosphereAnisotropy;

// Volumetric emission constants
uniform float volumeEmission;
uniform vec3 volumeEmissionColorRGB;

#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5
#define RADIANCE_EPSILON 1.0e-5

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2

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

float averageComponent(in vec3 v)
{
    return (v.r + v.g + v.b)/3.0;
}


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphereCosineWeighted(inout vec4 rnd, inout float pdf)
{
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}

#ifdef HAS_ATMOSPHERE
vec3 sampleSphereUniformly(inout vec4 rnd)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    // pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}
#endif

vec3 SURFACE_DIFFUSE_REFL_RGB(in vec3 X, in vec3 nW, in vec3 woW)
{
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, woW);
    return reflRGB;
}

#ifdef HAS_VOLUME_EMISSION
vec3 VOLUME_EMISSION_EVAL(in vec3 X)
{
    vec3 emission = VOLUME_EMISSION(volumeEmissionColorRGB * volumeEmission, X);
    return emission;
}
#endif // HAS_VOLUME_EMISSION

#ifdef HAS_SURFACE_NORMALMAP
vec3 perturbNormal(in vec3 X, in Basis basis, int material)
{
    if (material==MAT_SURFA)
    {
        vec3 nL = SURFACE_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL), basis);
    }
    return basis.nW;
}
#endif

//////////////////////////////////////////////
// Sky
//////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dirW)
{
    float rot_phi   = M_PI*envMapPhiRotation/180.0;
    float rot_theta = M_PI*envMapThetaRotation/180.0;
    vec3 sky_pole = vec3(sin(rot_theta), cos(rot_theta), 0.0);
    Basis sky_basis = makeBasis(sky_pole);
    vec3 dirL = worldToLocal(dirW, sky_basis);
    float phi = atan(dirL.y, dirL.x) + M_PI + rot_phi;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = acos(dirL.z);
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 RGB = vec3(1.0);
    if (haveEnvMap)
        RGB = texture(envMap, vec2(u,v)).rgb;
    float t = dot(dirW, sky_pole);
    float tt = envMapTransitionAngle/180.0;
    RGB *= skyPower * mix(skyTintDown, skyTintUp, smoothstep(-tt, tt, t));
    return RGB;
}

vec3 environmentRadiance(in vec3 dir)
{
    vec3 RGB_sky = environmentRadianceRGB(dir);
    return RGB_sky;
}

vec3 sampleSkyAtSurface(in Basis basis, inout vec4 rnd,
                        inout vec3 woutputW, inout float pdfDir)
{
    vec3 woutputL = sampleHemisphereCosineWeighted(rnd, pdfDir);
    woutputW = localToWorld(woutputL, basis);
    return environmentRadiance(woutputW);
}

#ifdef HAS_ATMOSPHERE
vec3 sampleSkyInVolume(inout vec4 rnd,
                       inout vec3 woutputW, inout float pdfDir)
{
    if (skyPower<RADIANCE_EPSILON)
        return vec3(0.0);
    woutputW = sampleSphereUniformly(rnd);
    pdfDir = 1.0/(4.0*M_PI);
    return environmentRadiance(woutputW);
}
#endif

//////////////////////////////////////////////
// Sun
//////////////////////////////////////////////

Basis sunBasis;

vec3 sampleSunDir(inout vec4 rnd, inout float pdfDir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    float theta = theta_max * sqrt(rand(rnd));
    float costheta = cos(theta);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    pdfDir = 1.0/solid_angle;
    return localToWorld(vec3(x, y, z), sunBasis);
}

float pdfSun(in vec3 dir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return 0.0;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    return 1.0/solid_angle;
}

vec3 sunRadiance(in vec3 dir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return vec3(0.0);
    vec3 RGB_sun = sunPower * sunColor;
    return RGB_sun;
}

vec3 sampleSunAtSurface(in Basis basis, inout vec4 rnd,
                        inout vec3 woutputW, inout float pdfDir)
{
    woutputW = sampleSunDir(rnd, pdfDir);
    vec3 woutputL = worldToLocal(woutputW, basis);
    if (woutputL.z < 0.0) return vec3(0.0);
    return sunRadiance(woutputW);
}

#ifdef HAS_ATMOSPHERE
vec3 sampleSunInVolume(inout vec4 rnd,
                       inout vec3 woutputW, inout float pdfDir)
{
    if (sunPower<RADIANCE_EPSILON)
        return vec3(0.0);
    woutputW = sampleSunDir(rnd, pdfDir);
    return sunRadiance(woutputW);
}
#endif


//////////////////////////////////////////////
// Atmosphere
//////////////////////////////////////////////

#ifdef HAS_ATMOSPHERE

vec3 VOLUME_SCATTERING_EVAL()
{
    return atmosphereScatteringCoeffRGB/(atmosphereMFP*lengthScale);
}

vec3 VOLUME_EXTINCTION_EVAL()
{
    return atmosphereExtinctionCoeffRGB/(atmosphereMFP*lengthScale);
}

#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }

bool atmosphereHit( in vec3 rayPos, in vec3 rayDir,
                    inout float t0, inout float t1 )
{
    vec3 dL = 1.0/rayDir;
    vec3 lo = (atmosphereBoundsMin*lengthScale - rayPos) * dL;
    vec3 hi = (atmosphereBoundsMax*lengthScale - rayPos) * dL;
    sort2(lo, hi);
    bool hit = !( lo.x>hi.y || lo.y>hi.x || lo.x>hi.z || lo.z>hi.x || lo.y>hi.z || lo.z>hi.y );
    t0 = max(max(lo.x, lo.y), lo.z);
    t1 = min(min(hi.x, hi.y), hi.z);
    return hit;
}

bool atmosphereSegment(in vec3 pW, in vec3 rayDir, float segmentLength, // input segment
                       inout float t0, inout float t1)                  // portion of input segment within atmosphere
{
    if ( !atmosphereHit(pW, rayDir, t0, t1) )
        return false;
    t0 = max(0.0, t0);
    if (t1<0.0 || segmentLength<t0)
        return false;
    t1 = min(segmentLength, t1);
    return true;
}

float phaseFunction(float mu, float anisotropy)
{
    float g = anisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}

#ifdef HAS_VOLUME_EMISSION
vec3 samplePhaseFunction(in vec3 rayDir, float anisotropy, inout vec4 rnd)
{
    float g = anisotropy;
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
#endif

vec3 directVolumeLighting(in vec3 pW, in vec3 rayDir, inout vec4 rnd)
{
    vec3 Ldirect = vec3(0.0);
    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        vec3 woutputW;
        float skyPdf;
        vec3 Li = sampleSkyInVolume(rnd, woutputW, skyPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            vec3 Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale);
            Li *= Tr;
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, skyPdf);
            }
        }
    }
    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        vec3 woutputW;
        float sunPdf;
        vec3 Li = sampleSunInVolume(rnd, woutputW, sunPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            vec3 Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale);
            Li *= Tr;
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, sunPdf);
            }
        }
    }
    return min(vec3(radianceClamp), Ldirect);
}

vec3 atmosphericInscatteringRadiance(in vec3 pW, in vec3 rayDir, in float segmentLength, inout vec4 rnd)
{
    vec3 Ls = vec3(0.0);
    vec3 extinction = VOLUME_EXTINCTION_EVAL();
    vec3 scattering = VOLUME_SCATTERING_EVAL();
    float extinction_norm = averageComponent(extinction);
    float scattering_norm = averageComponent(scattering);
    if ( extinction_norm < RADIANCE_EPSILON ||
         scattering_norm < RADIANCE_EPSILON ) return Ls;

    // Find sub-segment of supplied segment over which (homogeneous) atmosphere exists
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return Ls;

    // Sample a scattering point in [t0, t1] from a PDF proportional to the transmittance (normalized over [t0, t1])
    float T01 = exp(-(t1-t0)*extinction_norm);
    float t_scatter = t0 - log(1.0 - rand(rnd)*(1.0 - T01)) / extinction_norm; // sampled scatter distance
    vec3 pW_scatter = pW + t_scatter*rayDir;                                   // sampled scatter point
    vec3 opticalDepth_scatter = (t_scatter - t0) * extinction;         // optical depth from pW -> pW_scatter
    vec3 Tr_pW = exp(-opticalDepth_scatter);                           // transmittance from pW -> pW_scatter
    float invDistancePdf = (1.0 - T01) * exp(averageComponent(opticalDepth_scatter)) / extinction_norm; // PDF of scatter distance

    // Direct lighting
    vec3 Li = directVolumeLighting(pW_scatter, rayDir, rnd);
    Ls += Tr_pW * scattering * Li * invDistancePdf; // final estimator for scattered radiance

#ifdef HAS_VOLUME_EMISSION
    // Surface emission contribution
    vec3 woW = samplePhaseFunction(rayDir, atmosphereAnisotropy, rnd); // sample direction of scattered light (from phase function)
    vec3 pW_hit;
    int hitMaterial;
    bool hit = traceRay(pW_scatter, woW, pW_hit, hitMaterial, maxLengthScale);
    if (hit)
    {
        // (NB, phase function is the angle PDF, so the MC PDF denom. cancels it in the estimator)
        Ls += Tr_pW * scattering * VOLUME_EMISSION_EVAL(pW_hit) * invDistancePdf;
    }
#endif

    return Ls;
}

#endif // HAS_ATMOSPHERE


////////////////////////////////////////////////////////////////////////////////
// Ambient occlusion integrator
////////////////////////////////////////////////////////////////////////////////

// Return the amount of light transmitted (per-channel) along a known "free" segment (i.e. free of geometry)
vec3 transmittanceOverFreeSegment(in vec3 pW, in vec3 rayDir, float segmentLength)
{
    vec3 Tr = vec3(1.0);
#ifdef HAS_ATMOSPHERE
    vec3 extinction = VOLUME_EXTINCTION_EVAL();
    if (length(extinction) == 0.f)
        return vec3(1.0);
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return vec3(1.0);
    vec3 opticalDepth = (t1 - t0) * extinction;
    Tr = exp(-opticalDepth);
#endif
    return Tr;
}

// Return the amount of light transmitted (per-channel) along a given segment
// (zero if occluded by geometry).
vec3 transmittanceOverSegment(in vec3 pW, in vec3 rayDir, float segmentLength)
{
    vec3 pW_surface;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_surface, hitMaterial, segmentLength);
    if (hit)
        return vec3(0.0);
    else
        return transmittanceOverFreeSegment(pW, rayDir, segmentLength);
}

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
    vec4 rnd = texture(RngData, vTexCoord);
#ifdef INTERACTIVE_MODE
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
#endif

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    INIT();

    vec2 pixel = gl_FragCoord.xy;
    vec3 L = vec3(0.0);
    for (int n=0; n<__MAX_SAMPLES_PER_FRAME__; ++n)
    {
        // Jitter over pixel
        vec2 pixelj = pixel + (-0.5 + vec2(rand(rnd), rand(rnd)));
        vec3 primaryStart, primaryDir;
#ifdef HAS_CUSTOM_CAMERA
        CONSTRUCT_PRIMARY_RAY(pixelj, rnd, primaryStart, primaryDir);
#else
        constructPrimaryRay(pixelj, rnd, primaryStart, primaryDir);
#endif

        // Raycast to first hit point
        vec3 pW_hit;
        vec3 woW = -primaryDir;
        int hitMaterial;
        bool hit = traceRay(primaryStart, primaryDir, pW_hit, hitMaterial, maxLengthScale);

        // Closest-hit logic
        if (hit)
        {
            // Compute normal at hit point
            vec3 nW = normal(pW_hit, hitMaterial);
            Basis basis = makeBasis(nW);
#ifdef HAS_SURFACE_NORMALMAP
            nW = perturbNormal(pW, basis, hitMaterial);
            basis = makeBasis(nW);
#endif
            float rayLength = length(pW_hit - primaryStart);

            // Compute diffuse BSDF
            vec3 f = SURFACE_DIFFUSE_REFL_RGB(pW_hit, nW, woW) / M_PI;

            // Compute direct lighting
            vec3 Ldirect = vec3(0.0);
            vec3 dPw = 3.0*minLengthScale*nW;
            if (skyPower > RADIANCE_EPSILON)
            {
                // Sky
                float skyPdf;
                vec3 woutputW;
                vec3 Li = sampleSkyAtSurface(basis, rnd, woutputW, skyPdf);
                vec3 TrToLight = transmittanceOverSegment(pW_hit+dPw, woutputW, maxLengthScale);
                if (averageComponent(Li) > RADIANCE_EPSILON)
                    Ldirect += f * TrToLight * Li / max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, nW));
            }
            if (sunPower > RADIANCE_EPSILON)
            {
                // Sun
                float sunPdf;
                vec3 woutputW;
                vec3 Li = sampleSunAtSurface(basis, rnd, woutputW, sunPdf);
                vec3 TrToLight = transmittanceOverSegment(pW_hit+dPw, woutputW, maxLengthScale);
                if (averageComponent(Li) > RADIANCE_EPSILON)
                    Ldirect += f * TrToLight * Li / max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, nW));
            }
#ifdef HAS_ATMOSPHERE
            vec3 TrToHit = transmittanceOverFreeSegment(primaryStart, primaryDir, rayLength);
#else
            vec3 TrToHit = vec3(1.0);
#endif
            L += TrToHit * Ldirect;

#ifdef HAS_ATMOSPHERE
            // Add term for single-scattering in the homogeneous atmosphere over the segment up to the closest hit
            vec3 Lscatter = atmosphericInscatteringRadiance(primaryStart, primaryDir, rayLength, rnd);
            L += Lscatter;
#endif

#ifdef HAS_VOLUME_EMISSION
            // Add volumetric emission at the surface point, if present (treating it as an isotropic radiance field)
            L += TrToHit * VOLUME_EMISSION_EVAL(pW_hit);
#endif
        }

        // Miss logic
        else
        {
            vec3 TrToInfinity = vec3(1.0);
#ifdef HAS_ATMOSPHERE
            // Add term for single-inscattering in the homogeneous atmosphere over the segment to infinity
            vec3 Lscatter = atmosphericInscatteringRadiance(primaryStart, primaryDir, maxLengthScale, rnd);
            L += Lscatter;
            TrToInfinity = transmittanceOverFreeSegment(primaryStart, primaryDir, maxLengthScale);
#endif
            L += TrToInfinity * environmentRadiance(primaryDir);
            if (sunVisibleDirectly)
                L += TrToInfinity * sunRadiance(primaryDir);
        }

    } // sample n of __MAX_SAMPLES_PER_FRAME__

    L /= float(__MAX_SAMPLES_PER_FRAME__);

    vec3 XYZ = rgbToXyz(L);

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

'lighttracer-fragment-shader': `#version 300 es
precision highp float;

// @todo: samplers


// Spotlights
uniform int spotlightCount;
uniform vec3 spotlightPosition0;
uniform vec3 spotlightPosition1;
uniform vec3 spotlightPosition2;
uniform vec3 spotlightPosition3;
uniform vec3 spotlightDirection0;
uniform vec3 spotlightDirection1;
uniform vec3 spotlightDirection2;
uniform vec3 spotlightDirection3;
uniform float spotlightRadius0;
uniform float spotlightRadius1;
uniform float spotlightRadius2;
uniform float spotlightRadius3;
uniform float spotlightAngle0;
uniform float spotlightAngle1;
uniform float spotlightAngle2;
uniform float spotlightAngle3;
uniform vec3 spotlightRadiance0;
uniform vec3 spotlightRadiance1;
uniform vec3 spotlightRadiance2;
uniform vec3 spotlightRadiance3;

// Camera parameters
uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

// Length scales
uniform float lengthScale;
uniform float minLengthScale;
uniform float maxLengthScale;

// Surface material parameters
uniform float metalRoughness;
uniform vec3 metalSpecAlbedoRGB;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoRGB;
uniform vec3 surfaceDiffuseAlbedoRGB;
uniform vec3 surfaceSpecAlbedoRGB;
uniform float surfaceRoughness;
uniform float surfaceIor;

// Volumetric/atmosphere constants
uniform float volumeEmission;
uniform vec3 volumeEmissionColorRGB;
uniform float atmosphereExtinction;
uniform vec3 atmosphereScatteringColorRGB;
uniform vec3 atmosphereAbsorptionColorRGB;
uniform float atmosphereAnisotropy;



void constructSpotlightRay(in vec3 position, in vec3 direction, float radius, float angle, inout vec4 rnd,
                           inout vec3 primaryStart, inout vec3 primaryDir)
{
    // Make emission cross-section circular
    float rPos = radius*sqrt(rand(rnd));
    float phiPos = 2.0*M_PI*rand(rnd);
    vec3 X = vec3(1.0, 0.0, 0.0);
    vec3 Z = vec3(0.0, 0.0, 1.0);
    vec3 u = cross(Z, direction);
    if ( length(u) < 1.0e-3 )
        u = cross(X, direction);

    u = normalize(u);
    vec3 v = cross(direction, u);
    primaryStart = position + rPos*(u*cos(phiPos) + v*sin(phiPos));

    // Emit in a cone with the given spread
    float spreadAngle = 0.5*abs(angle)*M_PI/180.0;
    float rDir = min(tan(spreadAngle), 1.0e6) * sqrt(rand(rnd));
    float phiDir = 2.0*M_PI*rand(seed);
    primaryDir = normalize(direction + rDir*(u*cos(phiDir) + v*sin(phiDir)));
}



RadianceType lightPathContribution(float wavelength_nm, in vec3 rgb, inout vec4 rnd)
{
    // For each spotlight, trace a path from the light through the scene,
    // with appropriate vertex factors at each surface hit.
    // At each vertex, connect to the first camera-path hit point (cW), producing a final contribution
    // due to each light.

    // Also, on each light path segment, connect to the first camera segment [c0, c1]
    // to generate a volumetric scattering contribution (if an atmosphere is present).
    // This requires sampling the distances along each segment.
    // (Ideally we would importance sample these distances according to their joint PDF, but a good start
    // is to equi-angular sample the camera-segment based on the light position, and the light-segment based on the camera position).

    // Perform light-trace starting from spot-lights, if present:
    for (int n=0; n<spotlightCount*__MAX_SAMPLES_PER_FRAME__; ++n)
    {
        vec3 Lspot;
        vec3 primaryStart, primaryDir;
        switch (n)
        {
            case 0: { constructSpotlightRay(spotlightPosition0, spotlightDirection0, spotlightRadius0, spotlightAngle0, rnd, primaryStart, primaryDir); Lspot = spotlightRadiance0; break; }
            case 1: { constructSpotlightRay(spotlightPosition1, spotlightDirection1, spotlightRadius1, spotlightAngle1, rnd, primaryStart, primaryDir); Lspot = spotlightRadiance1; break; }
            case 2: { constructSpotlightRay(spotlightPosition2, spotlightDirection2, spotlightRadius2, spotlightAngle2, rnd, primaryStart, primaryDir); Lspot = spotlightRadiance2; break; }
            case 3: { constructSpotlightRay(spotlightPosition3, spotlightDirection3, spotlightRadius3, spotlightAngle3, rnd, primaryStart, primaryDir); Lspot = spotlightRadiance3; break; }
        }

        vec3 pW = primaryStart;
        vec3 rayDir = primaryDir; // (parallel to light direction)
    }
}


void main()
{
    // Read current photon vertex and direction

    // Calculate hit (or miss)

    // Compute contribution to caustic image (image plane hit and radiance)

    // Update photon vertex and direction



}
`,

'normals-fragment-shader': `#version 300 es
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

#ifdef HAS_SURFACE_NORMALMAP
vec3 perturbNormal(in vec3 X, in Basis basis, int material)
{
    if (material==MAT_SURFA)
    {
        vec3 nL_perturbed = SURFACE_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL_perturbed), basis);
    }
    return basis.nW;
}
#endif

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

    vec3 colorXYZ = vec3(0.0);
    for (int n=0; n<__MAX_SAMPLES_PER_FRAME__; ++n)
    {
        // Jitter over pixel
        vec2 pixelj = pixel + (-0.5 + vec2(rand(rnd), rand(rnd)));
        vec3 primaryStart, primaryDir;
        constructPrimaryRay(pixelj, rnd, primaryStart, primaryDir);

        // Raycast to first hit point
        vec3 pW;
        vec3 woW = -primaryDir;
        int hitMaterial;
        bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxLengthScale);
        if (hit)
        {
            // Compute normal at hit point
            vec3 nW = normal(pW, hitMaterial);
            vec3 ngW = nW; // geometric normal
            Basis basis = makeBasis(nW);
#ifdef HAS_SURFACE_NORMALMAP
            nW = perturbNormal(pW, basis, hitMaterial);
            // If the incident ray lies below the hemisphere of the perturbed shading normal,
            // which can occur due to normal mapping, apply the "Flipping hack" to prevent artifacts
            // (see Schüßler, "Microfacet-based Normal Mapping for Robust Monte Carlo Path Tracing")
            if ((dot(nW, primaryDir) > 0.0) && (hitMaterial != MAT_DIELE))
                nW = 2.0*ngW*dot(ngW, nW) - nW;
#endif
            colorXYZ += rgbToXyz(0.5*(nW+vec3(1.0)));
        }
        else
        {
            if (envMapVisible) colorXYZ += environmentRadianceXYZ(primaryDir);
            else               colorXYZ += vec3(0.0);
        }
    }
    colorXYZ /= float(__MAX_SAMPLES_PER_FRAME__);

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

'pathtracer-fragment-shader': `#version 300 es
precision highp float;

// Samplers
uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D iorTex;           // 4
uniform sampler2D kTex;             // 5 (for metals)
uniform sampler2D envMap;           // 6
in vec2 vTexCoord;

layout(location = 0) out vec4 gbuf_rad;
layout(location = 1) out vec4 gbuf_rng;

// Camera parameters
uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

// Pathtracing parameters
uniform float filterRadius;
uniform float radianceClamp;
uniform float skipProbability;
uniform float shadowStrength;
uniform bool maxStepsIsMiss;
uniform int wavelengthSamples;

// Length scales
uniform float lengthScale;
uniform float minLengthScale;
uniform float maxLengthScale;

// Sky parameters
uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapPhiRotation;
uniform float envMapThetaRotation;
uniform float envMapTransitionAngle;
uniform float skyPower;
uniform vec3 skyTintUp;
uniform vec3 skyTintDown;

// Sun parameters
uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;

// Sphere light
uniform vec3 sphereLightPosition;
uniform float sphereLightRadius;
uniform float sphereLightPower;
uniform vec3 sphereLightColor;

// Surface material parameters
uniform float metalRoughness;
uniform vec3 metalSpecAlbedoRGB;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoRGB;
uniform vec3 surfaceDiffuseAlbedoRGB;
uniform vec3 surfaceSpecAlbedoRGB;
uniform float surfaceRoughness;
uniform float surfaceIor;
uniform float subsurface;
uniform vec3 subsurfaceAlbedoRGB;
uniform float subsurfaceMFP;
uniform float subsurfaceAnisotropy;
uniform float subsurfaceDiffuseWeight;
uniform int maxSSSSteps;

// Atmosphere constants
uniform vec3 atmosphereBoundsMin;
uniform vec3 atmosphereBoundsMax;
uniform vec3 atmosphereExtinctionCoeffRGB;
uniform vec3 atmosphereScatteringCoeffRGB;
uniform float atmosphereMFP;
uniform float atmosphereAnisotropy;

// Fog
uniform bool fogEnable;
uniform float fogMaxOpticalDepth;
uniform float fogMFP;
uniform vec3 fogEmission;
uniform vec3 fogTint;

// Volumetric emission constants
uniform float volumeEmission;
uniform vec3 volumeEmissionColorRGB;


//////////////////////////////////////////////////////////////
// Defines
//////////////////////////////////////////////////////////////
#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5
#define RADIANCE_EPSILON 1.0e-6
#define M_PI 3.141592653589793

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2

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

#ifdef DISPERSION_ENABLED
#define RadianceType float
#else
#define RadianceType vec3
#endif

/////////////////////////////////////////////////////////////////////////
// Basis transforms
/////////////////////////////////////////////////////////////////////////

vec3 safe_normalize(vec3 N)
{
    float l = length(N);
    return N/max(l, DENOM_TOLERANCE);
}

float cosTheta2(in vec3 nLocal) { return nLocal.z*nLocal.z; }
float cosTheta(in vec3 nLocal)  { return nLocal.z; }
float sinTheta2(in vec3 nLocal) { return 1.0 - cosTheta2(nLocal); }
float sinTheta(in vec3 nLocal)  { return sqrt(max(0.0, sinTheta2(nLocal))); }
float tanTheta2(in vec3 nLocal) { float ct2 = cosTheta2(nLocal); return max(0.0, 1.0 - ct2) / max(ct2, 1.0e-7); }
float tanTheta(in vec3 nLocal)  { return sqrt(max(0.0, tanTheta2(nLocal))); }

struct Basis
{
    // right-handed coordinate system
    vec3 nW; // aligned with the z-axis in local space
    vec3 tW; // aligned with the x-axis in local space
    vec3 bW; // aligned with the y-axis in local space
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
    basis.tW = safe_normalize(basis.tW);
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

#ifdef DISPERSION_ENABLED
vec3 xyzToRgb(vec3 XYZ)
{
    // (Assuming RGB in sRGB color space)
    vec3 RGB;
    RGB.r =  3.2404542*XYZ.x - 1.5371385*XYZ.y - 0.4985314*XYZ.z;
    RGB.g = -0.9692660*XYZ.x + 1.8760108*XYZ.y + 0.0415560*XYZ.z;
    RGB.b =  0.0556434*XYZ.x - 0.2040259*XYZ.y + 1.0572252*XYZ.z;
    // deal with out-of-gamut RGB.
    float delta = -min(0.0, min(min(RGB.r, RGB.g), RGB.b));
    RGB.r += delta;
    RGB.g += delta;
    RGB.b += delta;
    // normalize
    float sum = RGB.r + RGB.g + RGB.b;
    RGB /= sum;
    return RGB;
}
#endif

vec3 rgbToXyz(in vec3 RGB)
{
    // (Assuming RGB in sRGB color space)
    vec3 XYZ;
    XYZ.x = 0.4124564*RGB.r + 0.3575761*RGB.g + 0.1804375*RGB.b;
    XYZ.y = 0.2126729*RGB.r + 0.7151522*RGB.g + 0.0721750*RGB.b;
    XYZ.z = 0.0193339*RGB.r + 0.1191920*RGB.g + 0.9503041*RGB.b;
    return XYZ;
}

#ifdef DISPERSION_ENABLED
// In dispersive rendering, takes the RGB of an albedo-like quantity (i.e. a desired color),
// and the rgb color matching functions at the wavelength of the current monochromatic beam,
// and returns the corresponding scalar albedo at this wavelength.
float rgbToAlbedo(in vec3 RGB, in vec3 rgb)
{
    return dot(RGB, rgb);
}
#else
vec3 rgbToAlbedo(in vec3 RGB, in vec3 rgb)
{
    return RGB;
}
#endif

float maxComponent(in vec3 v)
{
    return max(max(v.r, v.g), v.b);
}

float maxComponent(float f)
{
    return f;
}

float averageComponent(in vec3 v)
{
    return (v.r + v.g + v.b)/3.0;
}

float averageComponent(float f)
{
    return f;
}

int randomChannel(inout vec4 rnd)
{;
    float r = rand(rnd);
    if      (r<0.33333333) return 0;
    else if (r<0.66666666) return 1;
    else                   return 2;
}

float getChannel(in vec3 color, int channel)
{
    if      (channel==0) return color.r;
    else if (channel==1) return color.g;
    else                 return color.b;
}

void setChannel(inout vec3 color, int channel, float value)
{
    if      (channel==0) color.r = value;
    else if (channel==1) color.g = value;
    else                 color.b = value;
}


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

/// Uniformly sample a sphere
vec3 sampleSphereUniformly(inout vec4 rnd)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    // pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}

vec3 sampleHemisphereUniformly(inout vec4 rnd)
{
    float xi = rand(rnd);
    float r = sqrt(1.0 - xi*xi);
    float phi = 2.0 * M_PI * rand(rnd);
    float x = r * cos(phi);
    float y = r * sin(phi);
    float z = xi;
    // pdf = 1.0/(2.0*M_PI);
    return vec3(x, y, z);
}

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphereCosineWeighted(inout vec4 rnd, inout float pdf)
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
float pdfHemisphereCosineWeighted(in vec3 wiL)
{
    if (wiL.z < 0.0) return 0.0;
    return wiL.z / M_PI;
}

float powerHeuristic(const float a, const float b)
{
    return a/(a + b);
}


/////////////////////////////////////////////////////////////////////////
// Phase function
/////////////////////////////////////////////////////////////////////////

float phaseFunction(float mu, float anisotropy)
{
    float g = anisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}

vec3 samplePhaseFunction(in vec3 rayDir, float anisotropy, inout vec4 rnd)
{
    float g = anisotropy;
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
    return safe_normalize(vec3(sinThetaM*cosPhiM, sinThetaM*sinPhiM, cosThetaM));
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

// cosi         = magnitude of the cosine of the incident ray angle to the normal
// eta_ti       = ratio et/ei of the transmitted IOR (et) and incident IOR (ei)
float fresnelDielectricReflectance(in float cosi, in float eta_ti)
{
    float c = cosi;
    float g = eta_ti*eta_ti - 1.0 + c*c;
    if (g > 0.0)
    {
        g = sqrt(g);
        float A = (g-c) / (g+c);
        float B = (c*(g+c) - 1.0) / (c*(g-c) + 1.0);
        return 0.5*A*A*(1.0 + B*B);
    }
    return 1.0; // total internal reflection
}

/////////////////////////////////////////////////////////////////////////
// BRDF functions
/////////////////////////////////////////////////////////////////////////

#ifdef HAS_GEOMETRY

// ****************************        Dielectric        ****************************

#ifdef HAS_DIELECTRIC

// Given the direction (wt) of a light beam transmitted through a plane dielectric interface
// with the given normal (n) in any orientation, and the ratio eta_ti = et/ei of the transmitted IOR (et) and incident IOR (ei),
// compute the direction of the incident light beam (wi). Returns false if no such beam exists, due to total internal reflection.
bool refraction_given_transmitted_direction(in vec3 n, in float eta_ti, in vec3 wt, inout vec3 wi)
{
    float wtn = dot(wt, n);
    float disciminant = 1.0 - eta_ti*eta_ti*(1.0 - wtn*wtn);
    if (disciminant < 0.0) return false;
    wi = eta_ti*wt - n*sign(wtn)*(eta_ti*abs(wtn) - sqrt(disciminant));
    return true;
}

// Given the direction (wi) of a light beam incident on a plane dielectric interface
// with the given normal (n) in any orientation, and the ratio eta_ti = et/ei of the transmitted IOR (et) and incident IOR (ei),
// compute the direction of the transmitted light beam (wt). Returns false if no such beam exists, due to total internal reflection.
bool refraction_given_incident_direction(in vec3 n, in float eta_ti, in vec3 wi, inout vec3 wt)
{
    float win = dot(wi, n);
    float disciminant = 1.0 - (1.0 - win*win)/(eta_ti*eta_ti);
    if (disciminant < 0.0) return false;
    wt = wi/eta_ti - n*sign(win)*(abs(win)/eta_ti - sqrt(disciminant));
    return true;
}

RadianceType DIELECTRIC_SPEC_REFL_EVAL(in vec3 X, in vec3 winputL, in Basis basis, in vec3 rgb)
{
    vec3 winputW = localToWorld(winputL, basis);
    vec3 reflRGB = DIELECTRIC_SPECULAR_REFLECTANCE(dieleSpecAlbedoRGB, X, basis.nW, winputW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType DIELECTRIC_ABSORPTION_EVAL(in vec3 X, in vec3 rgb)
{
    vec3 absorptionRGB = DIELECTRIC_ABSORPTION(dieleAbsorptionRGB, X);
    return rgbToAlbedo(absorptionRGB, rgb);
}

RadianceType evaluateDielectric( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb, bool fromCamera )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(woutputL) * cosTheta(winputL) > 0.0;
    RadianceType dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, winputL, basis, rgb);
    vec3 beamIncidentL = fromCamera ? -woutputL : -winputL;
    vec3 beamOutgoingL = fromCamera ?  winputL  :  woutputL;
    bool entering = (beamIncidentL.z < 0.0);
    vec3 h;
    if (reflected)
    {
        // Compute reflection half-vector
        h = safe_normalize(woutputL + winputL);
    }
    else
    {
        // Compute refraction half-vector
        if (entering)
            h = safe_normalize(beamIncidentL - ior*beamOutgoingL);
        else
            h = safe_normalize(ior*beamIncidentL - beamOutgoingL);
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(-beamIncidentL, beamOutgoingL, h, roughness);
    RadianceType f;
    float eta_ti = entering ? ior : 1.0/ior; // eta_ti = et/ei
    RadianceType Fr = dielectricAlbedo * fresnelDielectricReflectance(abs(winputL.z), eta_ti);
    if (reflected)
    {
        f = dielectricAlbedo * Fr * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
    }
    else
    {
        float im = dot(-beamIncidentL, h);
        float om = dot(beamOutgoingL, h);
        float eta_ti = entering ? ior : 1.0/ior; // eta_it = ei/et
        float sqrtDenom = im + eta_ti*om;
        float dwh_dwo = eta_ti*eta_ti * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        f = (RadianceType(1.0) - Fr) * G * D * abs(im) * dwh_dwo / max(abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
    }
    return f;
}

float pdfDielectric( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb, bool fromCamera )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(woutputL) * cosTheta(winputL) > 0.0;
    RadianceType dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, winputL, basis, rgb);
    vec3 beamIncidentL = fromCamera ? -woutputL : -winputL;
    vec3 beamOutgoingL = fromCamera ?  winputL  :  woutputL;
    bool entering = (beamIncidentL.z < 0.0);
    float eta_ti = entering ? ior : 1.0/ior; // eta_it = ei/et
    RadianceType Fr = dielectricAlbedo * fresnelDielectricReflectance(abs(winputL.z), eta_ti);
    vec3 h;
    float dwh_dwo;
    float pdf;
    if (reflected)
    {
        h = safe_normalize(woutputL + winputL);
        dwh_dwo = 1.0 / max(4.0*abs(dot(winputL, h)), DENOM_TOLERANCE);
        pdf = averageComponent(Fr);
    }
    else
    {
        // Compute refraction half-vector
        if (entering)
            h = safe_normalize(beamIncidentL - ior*beamOutgoingL);
        else
            h = safe_normalize(ior*beamIncidentL - beamOutgoingL);
        if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
        float eta_ti = entering ? ior : 1.0/ior; // eta_it = ei/et
        float im = dot(-beamIncidentL, h);
        float om = dot(beamOutgoingL, h);
        float sqrtDenom = im + eta_ti*om;
        dwh_dwo = eta_ti*eta_ti * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        pdf = 1.0 - averageComponent(Fr);
    }
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    pdf *= microfacetPDF(h, roughness);
    return abs(pdf * dwh_dwo);
}

RadianceType sampleDielectric( in vec3 X, in Basis basis, in vec3 winputL, in float wavelength_nm, in vec3 rgb, bool fromCamera,
                               inout vec3 woutputL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_DIELE(wavelength_nm);
    RadianceType dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, winputL, basis, rgb);
    float eta_ti_refl = (winputL.z >= 0.0) ? ior : 1.0/ior; // et/ei on reflection of incident beam
    RadianceType Fr = dielectricAlbedo * fresnelDielectricReflectance(abs(winputL.z), eta_ti_refl);
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    float microPDF = microfacetPDF(m, roughness);
    float reflectProb = averageComponent(Fr);
    if (rand(rnd) < reflectProb)
    {
        // Compute specularly reflected ray direction
        woutputL = -winputL + 2.0*dot(winputL, m)*m; // Compute incident direction by reflecting winputL about m
        if (woutputL.z * winputL.z < 0.0) woutputL *= -1.0; // flip if reflected ray direction in wrong hemisphere
        float D = microfacetEval(m, roughness);
        float G = smithG2(winputL, woutputL, m, roughness); // Shadow-masking function
        RadianceType f = Fr * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, m)), DENOM_TOLERANCE);
        pdfOut = microPDF * reflectProb * dwh_dwo; // Return total BRDF and corresponding pdf
        return f;
     }
     else
     {
        vec3 beamIncidentL = fromCamera ? -woutputL : -winputL;
        vec3 beamOutgoingL = fromCamera ?  winputL  :  woutputL;
        bool entering = (beamIncidentL.z < 0.0);
        float eta_ti = entering ? ior : 1.0/ior; // et/ei on transmission of incident beam
        if (fromCamera)
        {
            if ( !refraction_given_transmitted_direction(m, eta_ti, beamOutgoingL, beamIncidentL) ) return RadianceType(0.0); // total internal reflection occurred
            woutputL = -beamIncidentL;
        }
        else
        {
            beamIncidentL = -winputL;
            if ( !refraction_given_incident_direction(m, eta_ti, beamIncidentL, beamOutgoingL) ) return RadianceType(0.0); // total internal reflection occurred
            woutputL = beamOutgoingL;
        }

        float im = dot(-beamIncidentL, m);
        RadianceType Frm = RadianceType(fresnelDielectricReflectance(abs(im), eta_ti));
        RadianceType Trm = RadianceType(1.0) - dielectricAlbedo * Frm;  // Fresnel transmittance
        // Evaluate microfacet distribution for the sampled half direction
        vec3 wh = m; // refraction half-vector = m
        float D = microfacetEval(wh, roughness);
        float G = smithG2(-beamIncidentL, beamOutgoingL, wh, roughness); // Shadow-masking function
        float dwh_dwo; // Jacobian of the half-direction mapping
        {
            float om = dot(beamOutgoingL, m);
            float sqrtDenom = im + eta_ti*om;
            dwh_dwo = eta_ti*eta_ti * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        }
        RadianceType f = RadianceType(abs(im) * dwh_dwo * Trm * G * D) / max(abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
        pdfOut = (1.0-reflectProb) * microPDF * abs(dwh_dwo);
        return f;
     }
}

#endif

// ****************************        Metal        ****************************

#ifdef HAS_METAL

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

RadianceType METAL_SPEC_REFL_EVAL(in vec3 X, in vec3 winputL, in Basis basis, in vec3 rgb)
{
    vec3 winputW = localToWorld(winputL, basis);
    vec3 reflRGB = METAL_SPECULAR_REFLECTANCE(metalSpecAlbedoRGB, X, basis.nW, winputW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType evaluateMetal( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb)
{
    if (winputL.z<0.0) return RadianceType(0.0);
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(winputL.z, ior, k);
    vec3 h = normalize(woutputL + winputL); // Compute the reflection half-vector
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(winputL, woutputL, h, roughness);
    RadianceType specAlbedo = METAL_SPEC_REFL_EVAL(X, winputL, basis, rgb);
    RadianceType f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
    return f;
}

float pdfMetal( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb)
{
    if (winputL.z<0.0) return PDF_EPSILON;
    float ior = IOR_DIELE(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    vec3 h = safe_normalize(woutputL + winputL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float pdf = microfacetPDF(h, roughness) * dwh_dwo;
    return pdf;
}

RadianceType sampleMetal( in vec3 X, in Basis basis, in vec3 winputL, in float wavelength_nm, in vec3 rgb,
                          inout vec3 woutputL, inout float pdfOut, inout vec4 rnd )
{
    if (winputL.z<0.0) return RadianceType(0.0);
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(winputL.z, ior, k);
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    woutputL = -winputL + 2.0*dot(winputL, m)*m; // Compute woutputL by reflecting winputL about m
    if (woutputL.z<DENOM_TOLERANCE) woutputL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
    float D = microfacetEval(m, roughness);
    float G = smithG2(winputL, woutputL, m, roughness); // Shadow-masking function
    RadianceType specAlbedo = METAL_SPEC_REFL_EVAL(X, winputL, basis, rgb);
    RadianceType f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
    float dwh_dwo; // Jacobian of the half-direction mapping
    dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, m)), DENOM_TOLERANCE);
    pdfOut = microfacetPDF(m, roughness) * dwh_dwo;
    return f;
}

#endif

// ****************************        Surface        ****************************

#ifdef HAS_SURFACE

RadianceType SURFACE_DIFFUSE_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 winputW, in vec3 rgb)
{
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, winputW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType SUBSURFACE_ALBEDO_EVAL(in vec3 X, in vec3 nW, in vec3 rgb)
{
    vec3 reflRGB = SUBSURFACE_ALBEDO(subsurfaceAlbedoRGB, X, nW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType SURFACE_SPEC_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 winputW, in vec3 rgb)
{
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(surfaceSpecAlbedoRGB, X, nW, winputW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType evaluateSurface(in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb)
{
    if (winputL.z<0.0) return RadianceType(0.0);
    vec3 winputW = localToWorld(winputL, basis);
    RadianceType diffuseAlbedo = (1.0 - subsurface)*SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, winputW, rgb);
    RadianceType    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, winputW, rgb);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float Fr = fresnelDielectricReflectance(woutputL.z, ior);
    vec3 h = normalize(woutputL + winputL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(winputL, woutputL, h, roughness);
    RadianceType f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
    float E = fresnelDielectricReflectance(winputL.z, ior);
    f += (1.0 - E) * diffuseAlbedo/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb)
{
    if (winputL.z<0.0) return PDF_EPSILON;
    vec3 winputW = localToWorld(winputL, basis);
    RadianceType diffuseAlbedo = (1.0 - subsurface)*SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, winputW, rgb);
    RadianceType    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, winputW, rgb);
    float ior = surfaceIor;
    float E = fresnelDielectricReflectance(abs(winputL.z), ior);
    float specWeight    = (E      )*averageComponent(specAlbedo);
    float diffuseWeight = (1.0 - E)*averageComponent(diffuseAlbedo);
    float weightSum = max(specWeight + diffuseWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float diffusePdf = pdfHemisphereCosineWeighted(woutputL);
    vec3 h = safe_normalize(woutputL + winputL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

RadianceType sampleSurface(in vec3 X, in Basis basis, in vec3 winputL, in float wavelength_nm, in vec3 rgb,
                           inout vec3 woutputL, inout float pdfOut, inout vec4 rnd)
{
    if (winputL.z<0.0) return RadianceType(0.0);
    vec3 winputW = localToWorld(winputL, basis);
    RadianceType diffuseAlbedo = (1.0 - subsurface)*SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, winputW, rgb);
    RadianceType    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, winputW, rgb);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float E = fresnelDielectricReflectance(abs(winputL.z), ior);
    float specWeight    = (E      )*averageComponent(specAlbedo);
    float diffuseWeight = (1.0 - E)*averageComponent(diffuseAlbedo);
    float weightSum = max(specWeight + diffuseWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        woutputL = sampleHemisphereCosineWeighted(rnd, pdfOut);
        pdfOut *= (1.0-specProb);
        return (1.0 - E) * diffuseAlbedo/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        woutputL = -winputL + 2.0*dot(winputL, m)*m; // Compute woutputL by reflecting winputL about m
        pdfOut = specProb;
        if (woutputL.z<DENOM_TOLERANCE) return RadianceType(0.0);
        float D = microfacetEval(m, roughness);
        float G = smithG2(winputL, woutputL, m, roughness); // Shadow-masking function
        float Fr = fresnelDielectricReflectance(abs(winputL.z), ior);
        RadianceType f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, m)), DENOM_TOLERANCE);
        pdfOut *= microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
}

#endif

// ****************************        BSDF common interface        ****************************

RadianceType evaluateBsdf( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in int material, in float wavelength_nm, in vec3 rgb, bool fromCamera,
                           inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    evaluateSurface(X, basis, winputL, woutputL, wavelength_nm, rgb); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      evaluateMetal(X, basis, winputL, woutputL, wavelength_nm, rgb); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return evaluateDielectric(X, basis, winputL, woutputL, wavelength_nm, rgb, fromCamera); }
#endif
}

RadianceType sampleBsdf( in vec3 X, in Basis basis, in vec3 winputL, in int material, in float wavelength_nm, in vec3 rgb, bool fromCamera,
                         inout vec3 woutputL, inout float pdfOut, inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    sampleSurface(X, basis, winputL, wavelength_nm, rgb, woutputL, pdfOut, rnd); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      sampleMetal(X, basis, winputL, wavelength_nm, rgb, woutputL, pdfOut, rnd); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return sampleDielectric(X, basis, winputL, wavelength_nm, rgb, fromCamera, woutputL, pdfOut, rnd); }
#endif
}

float pdfBsdf( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in int material, in float wavelength_nm, in vec3 rgb, bool fromCamera )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    pdfSurface(X, basis, winputL, woutputL, wavelength_nm, rgb); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      pdfMetal(X, basis, winputL, woutputL, wavelength_nm, rgb); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return pdfDielectric(X, basis, winputL, woutputL, wavelength_nm, rgb, fromCamera); }
#endif
}

#endif // HAS_GEOMETRY

RadianceType VOLUME_SCATTERING_EVAL(in vec3 rgb)
{
    return rgbToAlbedo(atmosphereScatteringCoeffRGB/(atmosphereMFP*lengthScale), rgb);
}

RadianceType VOLUME_EXTINCTION_EVAL(in vec3 rgb)
{
    return rgbToAlbedo(atmosphereExtinctionCoeffRGB/(atmosphereMFP*lengthScale), rgb);
}

#ifdef HAS_VOLUME_EMISSION

RadianceType VOLUME_EMISSION_EVAL(in vec3 X, in vec3 rgb)
{
    vec3 emission = VOLUME_EMISSION(volumeEmissionColorRGB * volumeEmission, X);
    return rgbToAlbedo(emission, rgb);
}

#endif // HAS_VOLUME_EMISSION

///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit along specified ray
bool traceRay(in vec3 start, in vec3 dir,
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_INVAL;
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

    float InitialSign = sign(sdf);
    float t = InitialSign * sdf; // (always take the first step along the ray direction)
    if (abs(t)>=maxMarchDist) return false;
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
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf;
        if (abs(t)>=maxMarchDist) return false;
    }
    return !maxStepsIsMiss;
}

vec3 normal(in vec3 pW, int material)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = 2.0*minLengthScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 Xp = pW+e.xyy; vec3 Xn = pW-e.xyy;
    vec3 Yp = pW+e.yxy; vec3 Yn = pW-e.yxy;
    vec3 Zp = pW+e.yyx; vec3 Zn = pW-e.yyx;
    vec3 N;
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { N = vec3(   SDF_SURFACE(Xp) -    SDF_SURFACE(Xn),    SDF_SURFACE(Yp) -    SDF_SURFACE(Yn),    SDF_SURFACE(Zp) -    SDF_SURFACE(Zn)); return safe_normalize(N); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { N = vec3(     SDF_METAL(Xp) -      SDF_METAL(Xn),      SDF_METAL(Yp) -      SDF_METAL(Yn),      SDF_METAL(Zp) -      SDF_METAL(Zn)); return safe_normalize(N); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { N = vec3(SDF_DIELECTRIC(Xp) - SDF_DIELECTRIC(Xn), SDF_DIELECTRIC(Yp) - SDF_DIELECTRIC(Yn), SDF_DIELECTRIC(Zp) - SDF_DIELECTRIC(Zn)); return safe_normalize(N); }
#endif
}

#ifdef HAS_NORMALMAP
vec3 perturbNormal(in vec3 X, in Basis basis, int material)
{
#ifdef HAS_SURFACE_NORMALMAP
    if (material==MAT_SURFA)
    {
        vec3 nL_perturbed = SURFACE_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL_perturbed), basis);
    }
#endif
#ifdef HAS_METAL_NORMALMAP
    if (material==MAT_METAL)
    {
        vec3 nL_perturbed = METAL_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL_perturbed), basis);
    }
#endif
#ifdef HAS_DIELECTRIC_NORMALMAP
    if (material==MAT_DIELE)
    {
        vec3 nL_perturbed = DIELECTRIC_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL_perturbed), basis);
    }
#endif
    return basis.nW;
}
#endif


#ifdef HAS_ATMOSPHERE

#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }

bool atmosphereHit( in vec3 rayPos, in vec3 rayDir,
                    inout float t0, inout float t1 )
{
    vec3 dL = 1.0/rayDir;
    vec3 lo = (atmosphereBoundsMin*lengthScale - rayPos) * dL;
    vec3 hi = (atmosphereBoundsMax*lengthScale - rayPos) * dL;
    sort2(lo, hi);
    bool hit = !( lo.x>hi.y || lo.y>hi.x || lo.x>hi.z || lo.z>hi.x || lo.y>hi.z || lo.z>hi.y );
    t0 = max(max(lo.x, lo.y), lo.z);
    t1 = min(min(hi.x, hi.y), hi.z);
    return hit;
}

bool atmosphereSegment(in vec3 pW, in vec3 rayDir, float segmentLength, // input segment
                       inout float t0, inout float t1)                  // portion of input segment within atmosphere
{
    if ( !atmosphereHit(pW, rayDir, t0, t1) )
        return false;
    t0 = max(0.0, t0);
    if (t1<0.0 || segmentLength<t0)
        return false;
    t1 = min(segmentLength, t1);
    return true;
}

#endif

// Return the amount of light transmitted (per-channel) along a known "free" segment (i.e. free of geometry)
RadianceType transmittanceOverFreeSegment(in vec3 pW, in vec3 rayDir, float segmentLength, in vec3 rgb)
{
    RadianceType Tr = RadianceType(1.0);
#ifdef HAS_ATMOSPHERE
    RadianceType extinction = VOLUME_EXTINCTION_EVAL(rgb);
    if (length(extinction) == 0.f)
        return RadianceType(1.0);
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return RadianceType(1.0);
    RadianceType opticalDepth = (t1 - t0) * extinction;
    Tr = exp(-opticalDepth);
#endif
#ifdef HAS_FOG
    if (fogEnable)
    {
        vec3 fogAbsorption = vec3(1.0) - fogTint;
        vec3 opticalDepthFog = fogMaxOpticalDepth * (vec3(1.0) - exp(-segmentLength*fogAbsorption/fogMFP));
        Tr *= exp(-rgbToAlbedo(opticalDepthFog, rgb));
    }
#endif
    return Tr;
}

// Return the amount of light transmitted (per-channel) along a given segment
// (zero if occluded by geometry).
RadianceType transmittanceOverSegment(in vec3 pW, in vec3 rayDir, float segmentLength, in vec3 rgb)
{
    vec3 pW_surface;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_surface, hitMaterial, segmentLength);
    if (hit)
        return RadianceType(0.0);
    else
        return transmittanceOverFreeSegment(pW, rayDir, segmentLength, rgb);
}



////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
// Sky
//////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dirW)
{
    float rot_phi   = M_PI*envMapPhiRotation/180.0;
    float rot_theta = M_PI*envMapThetaRotation/180.0;
    vec3 sky_pole = vec3(sin(rot_theta), cos(rot_theta), 0.0);
    Basis sky_basis = makeBasis(sky_pole);
    vec3 dirL = worldToLocal(dirW, sky_basis);
    float phi = atan(dirL.y, dirL.x) + M_PI + rot_phi;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = acos(dirL.z);
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 RGB = vec3(1.0);
    if (haveEnvMap)
        RGB = texture(envMap, vec2(u,v)).rgb;
    float t = dot(dirW, sky_pole);
    float tt = envMapTransitionAngle/180.0;
    RGB *= skyPower * mix(skyTintDown, skyTintUp, smoothstep(-tt, tt, t));
    return RGB;
}

RadianceType environmentRadiance(in vec3 dir, in vec3 rgb)
{
    vec3 RGB_sky = environmentRadianceRGB(dir);
    return rgbToAlbedo(RGB_sky, rgb);
}

float skyPowerEstimate()
{
    return 2.0*M_PI * skyPower * 0.5*(maxComponent(skyTintUp) + maxComponent(skyTintDown));
}

RadianceType sampleSkyAtSurface(Basis basis, in vec3 rgb, inout vec4 rnd,
                                inout vec3 woutputL, inout vec3 woutputW, inout float pdfDir)
{
    if (skyPower<RADIANCE_EPSILON)
        return RadianceType(0.0);
    woutputL = sampleHemisphereCosineWeighted(rnd, pdfDir);
    woutputW = localToWorld(woutputL, basis);
    return environmentRadiance(woutputW, rgb);
}

RadianceType sampleSkyInVolume(in vec3 rgb, inout vec4 rnd,
                               inout vec3 woutputW, inout float pdfDir)
{
    if (skyPower<RADIANCE_EPSILON)
        return RadianceType(0.0);
    woutputW = sampleSphereUniformly(rnd);
    pdfDir = 1.0/(4.0*M_PI);
    return environmentRadiance(woutputW, rgb);
}

//////////////////////////////////////////////
// Sun
//////////////////////////////////////////////

vec3 sampleSunDir(inout vec4 rnd, inout float pdfDir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    float theta = theta_max * sqrt(rand(rnd));
    float costheta = cos(theta);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    pdfDir = 1.0/solid_angle;
    return localToWorld(vec3(x, y, z), sunBasis);
}

float pdfSun(in vec3 dir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return 0.0;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    return 1.0/solid_angle;
}

RadianceType sunRadiance(in vec3 dir, in vec3 rgb)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return RadianceType(0.0);
    vec3 RGB_sun = sunPower * sunColor;
    return rgbToAlbedo(RGB_sun, rgb);
}

RadianceType sampleSunAtSurface(Basis basis, in vec3 rgb, inout vec4 rnd,
                                inout vec3 woutputL, inout vec3 woutputW, inout float pdfDir)
{
    if (sunPower<RADIANCE_EPSILON)
        return RadianceType(0.0);
    woutputW = sampleSunDir(rnd, pdfDir);
    woutputL = worldToLocal(woutputW, basis);
    if (woutputL.z < 0.0) return RadianceType(0.0);
    return sunRadiance(woutputW, rgb);
}

RadianceType sampleSunInVolume(in vec3 rgb, inout vec4 rnd,
                               inout vec3 woutputW, inout float pdfDir)
{
    if (sunPower<RADIANCE_EPSILON)
        return RadianceType(0.0);
    woutputW = sampleSunDir(rnd, pdfDir);
    return sunRadiance(woutputW, rgb);
}

//////////////////////////////////////////////
// Sphere light
//////////////////////////////////////////////

#ifdef HAS_SPHERE_LIGHT

vec3 sampleSphereLightPosition(inout vec4 rnd, in vec3 pW,
                               inout float pdfArea)
{
    // Sample a point uniformly on the sphere-light hemisphere facing the vertex pW
    vec3 cp = normalize(pW - sphereLightPosition);
    Basis hemisphereBasis = makeBasis(cp);
    vec3 vL = sampleHemisphereUniformly(rnd);
    vec3 vW = localToWorld(vL, hemisphereBasis);
    pdfArea = 1.0 / (2.0*M_PI*sphereLightRadius*sphereLightRadius); // @todo: precompute
    return sphereLightPosition + sphereLightRadius*vW;
}

RadianceType sphereLightRadiance(in vec3 pW, in vec3 dir, in vec3 rgb,
                                 inout vec3 pHit)
{
    // does ray (pW, dir) hit the sphere?
    vec3 m = pW - sphereLightPosition;
    float b = dot(m, dir);
    float c = dot(m, m) - sphereLightRadius*sphereLightRadius;
    if (c>0.0 && b>0.0)
        return RadianceType(0.0);
    float discr = b*b - c;
    if (discr < 0.0)
        return RadianceType(0.0);
    float t = -b - sqrt(discr);
    t = max(t, 0.0);
    pHit = pW + t*dir;
    vec3 RGB_spherelight = sphereLightPower * sphereLightColor;
    return rgbToAlbedo(RGB_spherelight, rgb);
}

RadianceType sampleSphereLight(in vec3 rgb, inout vec4 rnd, in vec3 pW,
                               inout vec3 pHit, inout float pdfDir)
{
    if (sphereLightPower<RADIANCE_EPSILON)
        return RadianceType(0.0);
    float pdfArea;
    pHit = sampleSphereLightPosition(rnd, pW, pdfArea);
    vec3 dHit = pHit - pW;
    float dHit2 = dot(dHit, dHit);
    dHit = normalize(dHit);
    vec3 nHit = normalize(pHit - sphereLightPosition);
    float costheta = abs(dot(nHit, dHit));
    pdfDir = pdfArea * dHit2 / max(DENOM_TOLERANCE, costheta); // convert area-measure PDF to solid-angle-measure
    vec3 RGB_spherelight = sphereLightPower * sphereLightColor;
    return rgbToAlbedo(RGB_spherelight, rgb);
}

#endif

//////////////////////////////////////////////
// Direct lighting routines
//////////////////////////////////////////////

// Estimate direct radiance at the given surface vertex
RadianceType directSurfaceLighting(in vec3 pW, Basis basis, in vec3 winputW, in int material,
                                   float wavelength_nm, in vec3 rgb, inout vec4 rnd,
                                   inout float skyPdf, inout float sunPdf, inout float sphPdf)
{
    vec3 winputL = worldToLocal(winputW, basis);
    bool fromCamera = true; // camera path
    RadianceType Ldirect = RadianceType(0.0);
    vec3 woutputL, woutputW;

    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        RadianceType Li = sampleSkyAtSurface(basis, rgb, rnd, woutputL, woutputW, skyPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            RadianceType Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = max(PDF_EPSILON, pdfBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera));
                RadianceType f = evaluateBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera, rnd);
                float misWeight = powerHeuristic(skyPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }

    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        RadianceType Li = sampleSunAtSurface(basis, rgb, rnd, woutputL, woutputW, sunPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            RadianceType Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = max(PDF_EPSILON, pdfBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera));
                RadianceType f = evaluateBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera, rnd);
                float misWeight = powerHeuristic(sunPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }

#ifdef HAS_SPHERE_LIGHT
    // Sphere light
    if ((sphereLightPower > RADIANCE_EPSILON) &&
        (dot(pW - sphereLightPosition, basis.nW) <= sphereLightRadius)) // (is sphere visible from the surface point?)
    {
        vec3 pHit;
        RadianceType Li = sampleSphereLight(rgb, rnd, pW, pHit, sphPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            float hitDist = max(0.0, length(pHit - pW) - 3.0*minLengthScale);
            woutputW = normalize(pHit - pW);
            woutputL = worldToLocal(woutputW, basis);
            RadianceType Tr = transmittanceOverSegment(pW, woutputW, hitDist, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = max(PDF_EPSILON, pdfBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera));
                RadianceType f = evaluateBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera, rnd);
                float misWeight = powerHeuristic(sphPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, sphPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }
#endif

    return min(RadianceType(radianceClamp), Ldirect);
}

#ifdef HAS_ATMOSPHERE

RadianceType directVolumeLighting(in vec3 pW, in vec3 rayDir, in vec3 rgb, inout vec4 rnd)
{
    bool fromCamera = true; // camera path
    RadianceType Ldirect = RadianceType(0.0);

    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        vec3 woutputW;
        float skyPdf;
        RadianceType Li = sampleSkyInVolume(rgb, rnd, woutputW, skyPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            RadianceType Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, skyPdf);
            }
        }
    }

    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        vec3 woutputW;
        float sunPdf;
        RadianceType Li = sampleSunInVolume(rgb, rnd, woutputW, sunPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            RadianceType Tr = transmittanceOverSegment(pW, woutputW, maxLengthScale, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, sunPdf);
            }
        }
    }

#ifdef HAS_SPHERE_LIGHT
    // Sphere light
    if (sphereLightPower > RADIANCE_EPSILON)
    {
        float sphPdf;
        vec3 pHit;
        RadianceType Li = sampleSphereLight(rgb, rnd, pW, pHit, sphPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            float hitDist = max(0.0, length(pHit - pW) - 3.0*minLengthScale);
            vec3 woutputW = normalize(pHit - pW);
            RadianceType Tr = transmittanceOverSegment(pW, woutputW, hitDist, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                float PF = phaseFunction(dot(-rayDir, -woutputW), atmosphereAnisotropy); // evaluate phase function for scattering -woutputW -> -rayDir
                Ldirect += PF * Li/max(PDF_EPSILON, sphPdf);
            }
        }
    }
#endif

    return min(RadianceType(radianceClamp), Ldirect);
}

RadianceType atmosphericInscatteringRadiance(in vec3 pW, in vec3 rayDir, in float segmentLength, in vec3 rgb, inout vec4 rnd)
{
    RadianceType Ls = RadianceType(0.0);

#ifdef HAS_FOG
    // Fog "fake" contribution
    if (fogEnable)
    {
        RadianceType TrFog = transmittanceOverFreeSegment(pW, rayDir, segmentLength, rgb);
        Ls = (RadianceType(1.0) - TrFog) * rgbToAlbedo(fogEmission, rgb);
    }
#endif

    RadianceType extinction = VOLUME_EXTINCTION_EVAL(rgb);
    RadianceType scattering = VOLUME_SCATTERING_EVAL(rgb);
    float extinction_norm = averageComponent(extinction);
    float scattering_norm = averageComponent(scattering);
    if ( extinction_norm < RADIANCE_EPSILON ||
         scattering_norm < RADIANCE_EPSILON ) return Ls;

    // Find sub-segment of supplied segment over which (homogeneous) atmosphere exists
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return Ls;

    // Sample a scattering point in [t0, t1] from a PDF proportional to the transmittance (normalized over [t0, t1])
    float T01 = exp(-(t1-t0)*extinction_norm);
    float t_scatter = t0 - log(1.0 - rand(rnd)*(1.0 - T01)) / extinction_norm; // sampled scatter distance
    vec3 pW_scatter = pW + t_scatter*rayDir;                                   // sampled scatter point
    RadianceType opticalDepth_scatter = (t_scatter - t0) * extinction;         // optical depth from pW -> pW_scatter
    RadianceType Tr_pW = exp(-opticalDepth_scatter);                           // transmittance from pW -> pW_scatter
    float invDistancePdf = (1.0 - T01) * exp(averageComponent(opticalDepth_scatter)) / extinction_norm; // PDF of scatter distance

    // Direct lighting
    RadianceType Li = directVolumeLighting(pW_scatter, rayDir, rgb, rnd);
    Ls += Tr_pW * scattering * Li * invDistancePdf; // final estimator for scattered radiance

#ifdef HAS_VOLUME_EMISSION
    // Surface emission contribution
    vec3 woW = samplePhaseFunction(rayDir, atmosphereAnisotropy, rnd); // sample direction of scattered light (from phase function)
    vec3 pW_hit;
    int hitMaterial;
    bool hit = traceRay(pW_scatter, woW, pW_hit, hitMaterial, maxLengthScale);
    if (hit)
    {
        // (NB, phase function is the angle PDF, so the MC PDF denom. cancels it in the estimator)
        Ls += Tr_pW * scattering * VOLUME_EMISSION_EVAL(pW_hit, rgb) * invDistancePdf;
    }
#endif

    return Ls;
}

RadianceType atmosphericScatterSample(in vec3 pW, in vec3 rayDir, float segmentLength, in vec3 rgb, inout vec4 rnd,
                                      inout vec3 pW_scatter, inout vec3 woutputW)
{
    RadianceType extinction = VOLUME_EXTINCTION_EVAL(rgb);
    RadianceType scattering = VOLUME_SCATTERING_EVAL(rgb);
    float extinction_norm = averageComponent(extinction);

    // Find sub-segment of supplied segment over which (homogeneous) atmosphere exists
    float t0, t1;
    bool hitAtmosphere = atmosphereSegment(pW, rayDir, segmentLength, t0, t1);
    if (!hitAtmosphere)
        return RadianceType(0.0);

    // Sample a scattering point in [t0, t1] from a PDF proportional to the transmittance (normalized over [t0, t1])
    float T01 = exp(-(t1-t0)*extinction_norm);
    float t_scatter = t0 - log(1.0 - rand(rnd)*(1.0 - T01)) / extinction_norm; // sampled scatter distance
    pW_scatter = pW + t_scatter*rayDir;                                        // sampled scatter point
    RadianceType opticalDepth_scatter = (t_scatter - t0) * extinction;         // optical depth from pW -> pW_scatter
    RadianceType Tr_pW = exp(-opticalDepth_scatter);                           // transmittance from pW -> pW_scatter
    float invDistancePdf = (1.0 - T01) * exp(averageComponent(opticalDepth_scatter)) / extinction_norm; // PDF of scatter distance
    woutputW = samplePhaseFunction(rayDir, atmosphereAnisotropy, rnd); // sample direction of scattered light (from phase function)

    return scattering * Tr_pW * invDistancePdf; // throughput for scattered ray
}

#endif

float atmosphericScatteringProbability(in RadianceType Tr, in vec3 rgb)
{
#ifdef HAS_ATMOSPHERE
    RadianceType extinction = VOLUME_EXTINCTION_EVAL(rgb);
    RadianceType scattering = VOLUME_SCATTERING_EVAL(rgb);
    float albedo_norm = averageComponent(scattering/extinction);
    float scatter_prob = albedo_norm * (1.0 - averageComponent(Tr));
    return scatter_prob;
#else
    return 0.0;
#endif
}


///////////////////////////////////////////////////////////////
// SSS
///////////////////////////////////////////////////////////////

#define HARDMAX_SSS_STEPS 2048
#define MIN_SSS_STEPS_BEFORE_RR 4

#ifdef HAS_SURFACE
bool randomwalk_SSS(in vec3 pW, in Basis basis, inout vec4 rnd,
                    in float MFP, in RadianceType subsurfaceAlbedo, in RadianceType diffuseAlbedoEntry,
                    inout RadianceType walk_throughput, inout vec3 pExit)
{
    // Returns true on successful walk to an exit point.
    // Otherwise false indicating termination without finding a valid exit point.
    vec3 pWalk = pW;
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    pWalk -= dPw; // Displace walk into interior to avoid issues with trace numerics
    float pdfDir;
    vec3 dirwalkL = sampleHemisphereCosineWeighted(rnd, pdfDir);
    vec3 dirwalkW = -localToWorld(dirwalkL, basis); // (negative sign to start walk in interior hemisphere)

    // We assume a diffuse Lambertian lobe at the entry point
    // (either pure white, or with the diffuse albedo of the entry point, or somewhere in between)
    RadianceType f = mix(RadianceType(1.0), diffuseAlbedoEntry, subsurfaceDiffuseWeight) / M_PI;
    RadianceType fOverPdf = min(RadianceType(radianceClamp), f/max(PDF_EPSILON, pdfDir));
    RadianceType surface_entry_throughput = fOverPdf * abs(dot(dirwalkW, basis.nW));
    walk_throughput = RadianceType(1.0); //surface_entry_throughput; // update walk throughput due to entry in medium

    for (int n=0; n<HARDMAX_SSS_STEPS; ++n)
    {
        if (n>maxSSSSteps)
            break;
        float walk_step = -log(rand(rnd)) * MFP * lengthScale;
        vec3 pHit;
        int hitMaterial;
        bool hit = traceRay(pWalk, dirwalkW, pHit, hitMaterial, walk_step);
        if (hit)
        {
            // walk left the surface
            pExit = pHit;
            return true;
        }
        // Point remains within the surface, continue walking.
        // First, make a Russian-roulette termination decision (after a minimum number of steps has been taken)
        float termination_prob = 0.0;
        if (n > MIN_SSS_STEPS_BEFORE_RR)
        {
            float continuation_prob = clamp(10.0*maxComponent(walk_throughput), 0.0, 1.0);
            float termination_prob = 1.0 - continuation_prob;
            if (rand(rnd) < termination_prob)
                break;
            walk_throughput /= continuation_prob; // update walk throughput due to RR continuation
        }
        dirwalkW = samplePhaseFunction(dirwalkW, subsurfaceAnisotropy, rnd);
        pWalk += walk_step*dirwalkW;
        walk_throughput *= subsurfaceAlbedo; // update walk throughput due to scattering in medium
    }
    return false;
}

// Estimate radiance at the SSS exit point
RadianceType SSS_exit_radiance(in vec3 pW, Basis basis,
                               in vec3 rgb, inout vec4 rnd,
                               in RadianceType diffuseAlbedoExit,
                               inout float skyPdf, inout float sunPdf, inout float sphPdf)
{
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    RadianceType Ldirect = RadianceType(0.0);
    vec3 woutputL, woutputW;

    // We assume a diffuse Lambertian lobe at the exit point
    // (either pure white, or with the diffuse albedo of the exit point, or somewhere in between)
    RadianceType f = mix(RadianceType(1.0), diffuseAlbedoExit, subsurfaceDiffuseWeight) / M_PI;

    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        RadianceType Li = sampleSkyAtSurface(basis, rgb, rnd, woutputL, woutputW, skyPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            RadianceType Tr = transmittanceOverSegment(pW+dPw, woutputW, maxLengthScale, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = pdfHemisphereCosineWeighted(woutputL);
                float misWeight = powerHeuristic(skyPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }
    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        RadianceType Li = sampleSunAtSurface(basis, rgb, rnd, woutputL, woutputW, sunPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            RadianceType Tr = transmittanceOverSegment(pW+dPw, woutputW, maxLengthScale, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = pdfHemisphereCosineWeighted(woutputL);
                float misWeight = powerHeuristic(sunPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }
#ifdef HAS_SPHERE_LIGHT
    // Sphere light
    if ((sphereLightPower > RADIANCE_EPSILON) &&
        (dot(pW - sphereLightPosition, basis.nW) <= sphereLightRadius)) // (is sphere visible from the surface point?)
    {
        vec3 pHit;
        RadianceType Li = sampleSphereLight(rgb, rnd, pW, pHit, sphPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            float hitDist = max(0.0, length(pHit - pW) - 3.0*minLengthScale);
            woutputW = normalize(pHit - pW);
            woutputL = worldToLocal(woutputW, basis);
            RadianceType Tr = transmittanceOverSegment(pW+dPw, woutputW, hitDist, rgb);
            Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-Tr));
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = pdfHemisphereCosineWeighted(woutputL);
                float misWeight = powerHeuristic(sphPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, sphPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }
#endif
    return min(RadianceType(radianceClamp), Ldirect);
}

#endif


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
    primaryDir = safe_normalize(camDir + s);
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
    primaryDir = safe_normalize(focalPlaneHit - lensPos);
}


RadianceType cameraPath(in vec3 primaryStart, in vec3 primaryDir, 
                        float wavelength_nm, in vec3 rgb, inout vec4 rnd)
{
    // Perform pathtrace starting from the camera lens to estimate the primary ray radiance, L
    RadianceType L = RadianceType(0.0);
    float misWeightSky = 1.0; // For MIS book-keeping
    float misWeightSun = 1.0; // For MIS book-keeping
#ifdef HAS_SPHERE_LIGHT
    float misWeightSph = 1.0; // For MIS book-keeping
#endif
    vec3 pW = primaryStart;
    vec3 rayDir = primaryDir; // (opposite to light beam direction)

#ifdef HAS_DIELECTRIC
    bool inDielectric = SDF_DIELECTRIC(primaryStart) < 0.0;
#endif
    RadianceType throughput = RadianceType(1.0);
    int atmosphere_scatters = 0;

    for (int vertex=0; vertex<=__MAX_BOUNCES__; ++vertex)
    {
        if (maxComponent(throughput) < THROUGHPUT_EPSILON) break;

        // Raycast along current propagation direction rayDir, from current vertex pW to pW_next
        vec3 pW_next;
        int hitMaterial;
        bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
        float rayLength = maxLengthScale;
        if (hit)
        {
            rayLength = length(pW_next - pW);
        }
        else
        {
            // This ray missed all geometry; add environment light term
            // (attenuated by transmittance through atmosphere to "infinity")
            // and terminate path
            float lightSegmentLength;
            RadianceType Tr;
#ifdef HAS_SPHERE_LIGHT
            // Check first for area-light (i.e. sphere light) hit
            vec3 pLightHit;
            RadianceType LiSph = (sphereLightPower<RADIANCE_EPSILON) ? RadianceType(0.0) : sphereLightRadiance(pW, rayDir, rgb, pLightHit);
            if (averageComponent(LiSph) > RADIANCE_EPSILON)
            {
                lightSegmentLength = length(pLightHit - pW);
                Tr = transmittanceOverFreeSegment(pW, rayDir, lightSegmentLength, rgb);
                float scatter_prob = atmosphericScatteringProbability(Tr, rgb);
                L += throughput * Tr * misWeightSph * LiSph / (1.0 - scatter_prob);
            }
            else
#endif
            // Add contribution from distant lights
            {
                lightSegmentLength = maxLengthScale;
                Tr = transmittanceOverFreeSegment(pW, rayDir, lightSegmentLength, rgb);
                float scatter_prob = atmosphericScatteringProbability(Tr, rgb);
                if (averageComponent(Tr) > RADIANCE_EPSILON)
                {
                    if (!(vertex==0 && !envMapVisible)      && misWeightSky>0.0) L += throughput * Tr * misWeightSky * environmentRadiance(rayDir, rgb) / (1.0 - scatter_prob);
                    if (!(vertex==0 && !sunVisibleDirectly) && misWeightSun>0.0) L += throughput * Tr * misWeightSun * sunRadiance(rayDir, rgb)         / (1.0 - scatter_prob);
                }
            }

#ifdef HAS_ATMOSPHERE
            // Add term for single-inscattering in the homogeneous atmosphere over the segment to the light hit
            if (vertex<2) // (restrict to first two segments only, for efficiency)
                L += throughput * atmosphericInscatteringRadiance(pW, rayDir, lightSegmentLength, rgb, rnd);

            // Scatter in atmosphere depending on transmittance over segment
            if (atmosphere_scatters+1 < __MAX_ATMOSPHERE_SCATTERS__)
            {
                float scatter_prob = atmosphericScatteringProbability(Tr, rgb);
                if (scatter_prob > rand(rnd))
                {
                    vec3 pW_scatter;
                    vec3 woutputW;
                    throughput *= atmosphericScatterSample(pW, rayDir, lightSegmentLength, rgb, rnd, pW_scatter, woutputW) / max(DENOM_TOLERANCE, scatter_prob);
                    pW = pW_scatter;
                    rayDir = woutputW;
                    atmosphere_scatters++;
                    continue;
                }
            }
#endif
            // Ray escapes to infinity
            break;
        }

        RadianceType Tr = RadianceType(1.0); // transmittance (about to be calculated) over segment to next hit
        bool hitSphereLight = false;

#ifdef HAS_SPHERE_LIGHT
        // Add possible contribution due to ray hitting the sphere light
        vec3 pLightHit;
        RadianceType LiSph = (sphereLightPower<RADIANCE_EPSILON) ? RadianceType(0.0) : sphereLightRadiance(pW, rayDir, rgb, pLightHit);
        if (averageComponent(LiSph) > RADIANCE_EPSILON)
        {
            float lightSegmentLength = length(pLightHit - pW);
            if (lightSegmentLength < rayLength)
            {
                hitSphereLight = true;
                Tr = transmittanceOverFreeSegment(pW, rayDir, lightSegmentLength, rgb);
                L += throughput * Tr * misWeightSph * LiSph;
#ifdef HAS_ATMOSPHERE
                // Add term for single-inscattering in the homogeneous atmosphere over the segment up to the sphere light hit
                if (vertex<2) // (restrict to first two segments only, for efficiency)
                    L += throughput * atmosphericInscatteringRadiance(pW, rayDir, lightSegmentLength, rgb, rnd);
#endif
            }
        }
#endif

        // If the current ray lies inside a dielectric, apply Beer's law for absorption
#ifdef HAS_DIELECTRIC
        if (inDielectric)
        {
            RadianceType absorption = DIELECTRIC_ABSORPTION_EVAL(pW, rgb);
            throughput *= exp(-rayLength*absorption);
        }
#endif

        // Add term for single-scattering in the homogeneous atmosphere over the segment up to the next hit
#ifdef HAS_ATMOSPHERE
        if (vertex<2) // (restrict to first two segments only, for efficiency)
            L += throughput * atmosphericInscatteringRadiance(pW, rayDir, rayLength, rgb, rnd);

        if (!hitSphereLight)
            Tr = transmittanceOverFreeSegment(pW, rayDir, rayLength, rgb);

        // Scatter in atmosphere depending on transmittance over segment to next hit
        if (atmosphere_scatters+1 < __MAX_ATMOSPHERE_SCATTERS__)
        {
            float scatter_prob = atmosphericScatteringProbability(Tr, rgb);
            if (scatter_prob > rand(rnd))
            {
                vec3 pW_scatter;
                vec3 woutputW;
                throughput *= atmosphericScatterSample(pW, rayDir, rayLength, rgb, rnd, pW_scatter, woutputW) / max(DENOM_TOLERANCE, scatter_prob);
                pW = pW_scatter;
                rayDir = woutputW;
                atmosphere_scatters++;
                continue;
            }
            else
                throughput /= max(DENOM_TOLERANCE, 1.0 - scatter_prob);
        }

        if (hitSphereLight)
            break; // terminate on hitting light (as we assume the sphere light has zero albedo)

        // Attenuate throughput to next hit due to atmospheric attenuation (Beer's law)
        throughput *= transmittanceOverFreeSegment(pW, rayDir, rayLength, rgb);
#endif

        // This ray hit some geometry, so deal with the surface interaction.
        // First, compute the normal and thus the local vertex basis:
        pW = pW_next;
        vec3 nW = normal(pW, hitMaterial);
        vec3 ngW = nW; // geometric normal, needed for robust ray continuation
        Basis basis = makeBasis(nW);
#ifdef HAS_NORMALMAP
        nW = perturbNormal(pW, basis, hitMaterial);
        // If the incident ray lies below the hemisphere of the perturbed shading normal,
        // which can occur due to normal mapping, apply the "Flipping hack" to prevent artifacts
        // (see Schüßler, "Microfacet-based Normal Mapping for Robust Monte Carlo Path Tracing")
        if ((dot(nW, rayDir) > 0.0) && (hitMaterial != MAT_DIELE))
            nW = 2.0*ngW*dot(ngW, nW) - nW;
        basis = makeBasis(nW);
#endif

        // Make a binary choice whether to scatter at the surface, or do a subsurface random walk:
        bool do_subsurface_walk = false;
        float prob_sss = 0.0;
        RadianceType subsurfaceAlbedo;
#ifdef HAS_SURFACE
        if (subsurfaceMFP > 0.0)         // a) the surface must have non-zero subsurface MFP, and
        {
            subsurfaceAlbedo = SUBSURFACE_ALBEDO_EVAL(pW, basis.nW, rgb);
            float diffuse_weight    = (1.0 - subsurface)*averageComponent(SURFACE_DIFFUSE_REFL_EVAL(pW, basis.nW, -rayDir, rgb));
            float    spec_weight    = averageComponent(SURFACE_SPEC_REFL_EVAL(pW, basis.nW, -rayDir, rgb));
            float subsurface_weight = subsurface * averageComponent(subsurfaceAlbedo) * 3.0; // weight more highly due to higher variance
            float total_weight = diffuse_weight + spec_weight + subsurface_weight;
            prob_sss = subsurface_weight / (total_weight + DENOM_TOLERANCE);
            prob_sss = clamp(prob_sss, 0.0, 1.0-PDF_EPSILON);
            do_subsurface_walk = (rand(rnd) < prob_sss);
        }
#endif

        // Regular surface bounce case
        if (!do_subsurface_walk)
        {
            // Sample BSDF for the next bounce direction
            vec3 winputW = -rayDir; // winputW, points *towards* the incident direction
            vec3 winputL = worldToLocal(winputW, basis);
            vec3 woutputL; // woutputL, points *towards* the outgoing direction
            bool fromCamera = true; // camera path
            float bsdfPdf;
            RadianceType f = sampleBsdf(pW, basis, winputL, hitMaterial, wavelength_nm, rgb, fromCamera, woutputL, bsdfPdf, rnd);
            vec3 woutputW = localToWorld(woutputL, basis);

            // Detect dielectric transmission
#ifdef HAS_DIELECTRIC
            if (hitMaterial==MAT_DIELE && dot(winputW, ngW)*dot(woutputW, ngW) < 0.0)
                inDielectric = !inDielectric;
#endif
#ifdef HAS_VOLUME_EMISSION
            // Add volumetric emission at the surface point, if present (treating it as an isotropic radiance field)
            L += throughput * VOLUME_EMISSION_EVAL(pW, rgb);
#endif

            // Update ray direction to the BSDF-sampled direction
            rayDir = woutputW;

            // Prepare for tracing the direct lighting and continuation rays
            pW += ngW * sign(dot(rayDir, ngW)) * 3.0*minLengthScale; // perturb vertex into half-space of scattered ray

            // Add direct lighting term at current surface vertex
            float skyPdf = 0.0;
            float sunPdf = 0.0;
            float sphPdf = 0.0;
#ifdef HAS_DIELECTRIC
            if (!inDielectric)
#endif
                L += throughput * directSurfaceLighting(pW, basis, winputW, hitMaterial, wavelength_nm, rgb, rnd,
                                                        skyPdf, sunPdf, sphPdf);

            // Update path continuation throughput
            RadianceType fOverPdf = min(RadianceType(radianceClamp), f/max(PDF_EPSILON, bsdfPdf));
            RadianceType surface_throughput = fOverPdf * abs(dot(woutputW, nW)) / max(PDF_EPSILON, 1.0 - prob_sss);
            throughput *= surface_throughput;

            misWeightSky = powerHeuristic(bsdfPdf, skyPdf); // compute sky MIS weight for bounce ray
            misWeightSun = powerHeuristic(bsdfPdf, sunPdf); // compute sun MIS weight for bounce ray
#ifdef HAS_SPHERE_LIGHT
            misWeightSph = powerHeuristic(bsdfPdf, sphPdf); // compute sphere-light MIS weight for bounce ray
#endif
        }

#ifdef HAS_SURFACE
        // Do subsurface random walk to a new vertex
        else
        {
            RadianceType diffuseAlbedoEntry = SURFACE_DIFFUSE_REFL_EVAL(pW, basis.nW, -rayDir, rgb);
            RadianceType walk_throughput;
            vec3 pExit;
            basis.nW = ngW; // (use basis with geometric normal for SSS entry, to avoid surface acne)
            bool success = randomwalk_SSS(pW, basis, rnd, subsurfaceMFP, subsurfaceAlbedo, diffuseAlbedoEntry, walk_throughput, pExit);
            if (!success)
                break;

            // Update path vertex to exit point
            pW = pExit;

            // Compute updated normal and basis at exit point
            nW = normal(pExit, MAT_SURFA);
            Basis basis_exit = makeBasis(nW);

            // Sample direction of exit ray (assuming a diffuse Lambertian lobe with diffuse albedo of exit point)
            float bsdfPdf;
            vec3 dirExitL = sampleHemisphereCosineWeighted(rnd, bsdfPdf);
            vec3 woutputW = localToWorld(dirExitL, basis_exit);

            // Update ray direction to the SSS exit direction
            rayDir = woutputW;

            // Update throughput due to random walk
            throughput *= walk_throughput;

            // Add direct lighting term at exit vertex (assumed to be a diffuse lobe)
            float skyPdf = 0.0;
            float sunPdf = 0.0;
            float sphPdf = 0.0;
            RadianceType diffuseAlbedoExit = SURFACE_DIFFUSE_REFL_EVAL(pW, nW, -rayDir, rgb);
#ifdef HAS_DIELECTRIC
            if (!inDielectric)
#endif
                L += throughput * SSS_exit_radiance(pExit, basis_exit, rgb, rnd, diffuseAlbedoExit,
                                                    skyPdf, sunPdf, sphPdf);
            misWeightSky = powerHeuristic(bsdfPdf, skyPdf); // compute sky MIS weight for bounce ray
            misWeightSun = powerHeuristic(bsdfPdf, sunPdf); // compute sun MIS weight for bounce ray
#ifdef HAS_SPHERE_LIGHT
            misWeightSph = powerHeuristic(bsdfPdf, sphPdf); // compute sphere-light MIS weight for bounce ray
#endif

            // Update path continuation throughput
            RadianceType f = mix(RadianceType(1.0), diffuseAlbedoExit, subsurfaceDiffuseWeight) / M_PI;
            RadianceType fOverPdf = min(RadianceType(radianceClamp), f/max(PDF_EPSILON, bsdfPdf));
            RadianceType surface_exit_throughput = fOverPdf * abs(dot(woutputW, nW));
            throughput *= surface_exit_throughput / max(PDF_EPSILON, prob_sss);

            // Prepare for tracing the continuation ray
            pW += nW * sign(dot(rayDir, nW)) * 3.0*minLengthScale; // perturb vertex into half-space of scattered ray
        }
#endif

        // Bail out now if the path continuation throughput is below threshold
        if (maxComponent(throughput) < THROUGHPUT_EPSILON)
            break;
    }
    return L;
}

float sample_jitter(float xi)
{
    // sample from triangle filter, returning jitter in [-1, 1]
    return xi < 0.5 ? sqrt(2.0*xi) - 1.0 : 1.0 - sqrt(2.0 - 2.0*xi);
}

void pathtrace(vec2 pixel, vec4 rnd) // the current pixel
{
#ifdef DISPERSION_ENABLED
    // Sample photon wavelength via the inverse CDF of the emission spectrum.
    // We limit the sampled wavelengths to between about 420nm and 700nm to avoid color mismatch at low sample count
    float xi = 0.1+0.8*(0.5+floor(float(wavelengthSamples)*rand(rnd))) / float(wavelengthSamples);
    float w = texture(ICDF, vec2(xi, 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;
    vec3 xyz = texture(WavelengthToXYZ, vec2(w, 0.5)).rgb; // xyz CIE color matching functions
    vec3 rgb = xyzToRgb(xyz); // corresponding normalized rgb color matching functions
#else
    float wavelength_nm = 390.0 + (750.0 - 390.0)*0.5; // non-dispersive rendering uses mid-range wavelength
    vec3 rgb = vec3(0.0);
#endif

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    // Sample radiance of primary ray
    RadianceType L = RadianceType(0.0);
    for (int n=0; n<__MAX_SAMPLES_PER_FRAME__; ++n)
    {
        // Apply FIS to obtain pixel jitter about center in pixel units
        float jx = 0.5 * filterRadius * sample_jitter(rand(rnd));
        float jy = 0.5 * filterRadius * sample_jitter(rand(rnd));
        vec2 pixelj = pixel + vec2(jx, jy);

        // Compute world ray direction for this fragment
        vec3 primaryStart, primaryDir;
#ifdef HAS_CUSTOM_CAMERA
        CONSTRUCT_PRIMARY_RAY(pixelj, rnd, primaryStart, primaryDir);
#else
        constructPrimaryRay(pixelj, rnd, primaryStart, primaryDir);
#endif

        // Perform pathtrace to estimate the primary ray radiance, L
        L += cameraPath(primaryStart, primaryDir, wavelength_nm, rgb, rnd);
    }
    L /= float(__MAX_SAMPLES_PER_FRAME__);

    // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;

    // Compute tristimulus contribution from estimated radiance
#ifdef DISPERSION_ENABLED
    vec3 colorXYZ = xyz * L;
#else
    vec3 colorXYZ = rgbToXyz(L);
#endif

    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;
    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}

void main()
{
    vec4 rnd = texture(RngData, vTexCoord);
#ifdef INTERACTIVE_MODE
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
#endif

    INIT();
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

////////////////////////////////////////////////////////////////////////////////
// Picking program: computes first hit distance based on mouse location
////////////////////////////////////////////////////////////////////////////////

void constructPickRay(in vec2 pixel,
                      inout vec3 primaryStart, inout vec3 primaryDir)
{
    // Compute world ray direction for given fragment
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

'simplepathtracer-fragment-shader': `#version 300 es
precision highp float;

// Samplers
uniform sampler2D Radiance;         // 0 (IO buffer)
uniform sampler2D RngData;          // 1 (IO buffer)
uniform sampler2D WavelengthToXYZ;  // 2
uniform sampler2D ICDF;             // 3
uniform sampler2D iorTex;           // 4
uniform sampler2D kTex;             // 5 (for metals)
uniform sampler2D envMap;           // 6
in vec2 vTexCoord;

layout(location = 0) out vec4 gbuf_rad;
layout(location = 1) out vec4 gbuf_rng;

// Camera parameters
uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

// Pathtracing parameters
uniform float filterRadius;
uniform float radianceClamp;
uniform float skipProbability;
uniform bool maxStepsIsMiss;

// Length scales
uniform float lengthScale;
uniform float minLengthScale;
uniform float maxLengthScale;

// Sky parameters
uniform bool haveEnvMap;
uniform bool envMapVisible;
uniform float envMapPhiRotation;
uniform float envMapThetaRotation;
uniform float envMapTransitionAngle;
uniform float skyPower;
uniform vec3 skyTintUp;
uniform vec3 skyTintDown;

// Sun parameters
uniform float sunPower;
uniform float sunAngularSize;
uniform float sunLatitude;
uniform float sunLongitude;
uniform vec3 sunColor;
uniform vec3 sunDir;
uniform bool sunVisibleDirectly;

// Surface material parameters
uniform vec3 surfaceDiffuseAlbedoRGB;
uniform vec3 surfaceSpecAlbedoRGB;
uniform float surfaceRoughness;
uniform float surfaceIor;
uniform float subsurface;
uniform vec3 subsurfaceAlbedoRGB;
uniform float subsurfaceMFP;
uniform float subsurfaceAnisotropy;
uniform float subsurfaceDiffuseWeight;
uniform int maxSSSSteps;

// Volumetric emission constants
uniform float volumeEmission;
uniform vec3 volumeEmissionColorRGB;

// Fog
uniform bool fogEnable;
uniform float fogMaxOpticalDepth;
uniform float fogMFP;
uniform vec3 fogEmission;
uniform vec3 fogTint;


//////////////////////////////////////////////////////////////
// Defines
//////////////////////////////////////////////////////////////
#define DENOM_TOLERANCE 1.0e-7
#define PDF_EPSILON 1.0e-6
#define THROUGHPUT_EPSILON 1.0e-5
#define RADIANCE_EPSILON 1.0e-6
#define M_PI 3.141592653589793

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2

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

vec3 safe_normalize(vec3 N)
{
    float l = length(N);
    return N/max(l, DENOM_TOLERANCE);
}

float cosTheta2(in vec3 nLocal) { return nLocal.z*nLocal.z; }
float cosTheta(in vec3 nLocal)  { return nLocal.z; }
float sinTheta2(in vec3 nLocal) { return 1.0 - cosTheta2(nLocal); }
float sinTheta(in vec3 nLocal)  { return sqrt(max(0.0, sinTheta2(nLocal))); }
float tanTheta2(in vec3 nLocal) { float ct2 = cosTheta2(nLocal); return max(0.0, 1.0 - ct2) / max(ct2, 1.0e-7); }
float tanTheta(in vec3 nLocal)  { return sqrt(max(0.0, tanTheta2(nLocal))); }

struct Basis
{
    // right-handed coordinate system
    vec3 nW; // aligned with the z-axis in local space
    vec3 tW; // aligned with the x-axis in local space
    vec3 bW; // aligned with the y-axis in local space
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
    basis.tW = safe_normalize(basis.tW);
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

float maxComponent(in vec3 v)
{
    return max(max(v.r, v.g), v.b);
}

float maxComponent(float f)
{
    return f;
}

float averageComponent(in vec3 v)
{
    return (v.r + v.g + v.b)/3.0;
}

float averageComponent(float f)
{
    return f;
}


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

/// Uniformly sample a sphere
vec3 sampleSphereUniformly(inout vec4 rnd)
{
    float z = 1.0 - 2.0*rand(rnd);
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0*M_PI*rand(rnd);
    float x = cos(phi);
    float y = sin(phi);
    // pdf = 1.0/(4.0*M_PI);
    return vec3(x, y, z);
}

vec3 sampleHemisphereUniformly(inout vec4 rnd)
{
    float xi = rand(rnd);
    float r = sqrt(1.0 - xi*xi);
    float phi = 2.0 * M_PI * rand(rnd);
    float x = r * cos(phi);
    float y = r * sin(phi);
    float z = xi;
    // pdf = 1.0/(2.0*M_PI);
    return vec3(x, y, z);
}

// Do cosine-weighted sampling of hemisphere
vec3 sampleHemisphereCosineWeighted(inout vec4 rnd, inout float pdf)
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
float pdfHemisphereCosineWeighted(in vec3 wiL)
{
    if (wiL.z < 0.0) return 0.0;
    return wiL.z / M_PI;
}

float powerHeuristic(const float a, const float b)
{
    return a/(a + b);
}


/////////////////////////////////////////////////////////////////////////
// Phase function
/////////////////////////////////////////////////////////////////////////

float phaseFunction(float mu, float anisotropy)
{
    float g = anisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}

vec3 samplePhaseFunction(in vec3 rayDir, float anisotropy, inout vec4 rnd)
{
    float g = anisotropy;
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
    return safe_normalize(vec3(sinThetaM*cosPhiM, sinThetaM*sinPhiM, cosThetaM));
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

// ****************************        Surface        ****************************

#ifdef HAS_SURFACE

// cosi         = magnitude of the cosine of the incident ray angle to the normal
// eta_ti       = ratio et/ei of the transmitted IOR (et) and incident IOR (ei)
float fresnelDielectricReflectance(in float cosi, in float eta_ti)
{
    float c = cosi;
    float g = eta_ti*eta_ti - 1.0 + c*c;
    if (g > 0.0)
    {
        g = sqrt(g);
        float A = (g-c) / (g+c);
        float B = (c*(g+c) - 1.0) / (c*(g-c) + 1.0);
        return 0.5*A*A*(1.0 + B*B);
    }
    return 1.0; // total internal reflection
}

vec3 SURFACE_DIFFUSE_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 winputW, in vec3 rgb)
{
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, winputW);
    return reflRGB;
}

vec3 SUBSURFACE_ALBEDO_EVAL(in vec3 X, in vec3 nW, in vec3 rgb)
{
    vec3 reflRGB = SUBSURFACE_ALBEDO(subsurfaceAlbedoRGB, X, nW);
    return reflRGB;
}

vec3 SURFACE_SPEC_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 winputW, in vec3 rgb)
{
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(surfaceSpecAlbedoRGB, X, nW, winputW);
    return reflRGB;
}

vec3 evaluateSurface(in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb)
{
    if (winputL.z<0.0) return vec3(0.0);
    vec3 winputW = localToWorld(winputL, basis);
    vec3 diffuseAlbedo = (1.0 - subsurface)*SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, winputW, rgb);
    vec3    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, winputW, rgb);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float Fr = fresnelDielectricReflectance(woutputL.z, ior);
    vec3 h = normalize(woutputL + winputL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(winputL, woutputL, h, roughness);
    vec3 f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
    float E = fresnelDielectricReflectance(winputL.z, ior);
    f += (1.0 - E) * diffuseAlbedo/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in float wavelength_nm, in vec3 rgb)
{
    if (winputL.z<0.0) return PDF_EPSILON;
    vec3 winputW = localToWorld(winputL, basis);
    vec3 diffuseAlbedo = (1.0 - subsurface)*SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, winputW, rgb);
    vec3    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, winputW, rgb);
    float ior = surfaceIor;
    float E = fresnelDielectricReflectance(abs(winputL.z), ior);
    float specWeight    = (E      )*averageComponent(specAlbedo);
    float diffuseWeight = (1.0 - E)*averageComponent(diffuseAlbedo);
    float weightSum = max(specWeight + diffuseWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float diffusePdf = pdfHemisphereCosineWeighted(woutputL);
    vec3 h = safe_normalize(woutputL + winputL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

vec3 sampleSurface(in vec3 X, in Basis basis, in vec3 winputL, in float wavelength_nm, in vec3 rgb,
                           inout vec3 woutputL, inout float pdfOut, inout vec4 rnd)
{
    if (winputL.z<0.0) return vec3(0.0);
    vec3 winputW = localToWorld(winputL, basis);
    vec3 diffuseAlbedo = (1.0 - subsurface)*SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, winputW, rgb);
    vec3    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, winputW, rgb);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float E = fresnelDielectricReflectance(abs(winputL.z), ior);
    float specWeight    = (E      )*averageComponent(specAlbedo);
    float diffuseWeight = (1.0 - E)*averageComponent(diffuseAlbedo);
    float weightSum = max(specWeight + diffuseWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        woutputL = sampleHemisphereCosineWeighted(rnd, pdfOut);
        pdfOut *= (1.0-specProb);
        return (1.0 - E) * diffuseAlbedo/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        woutputL = -winputL + 2.0*dot(winputL, m)*m; // Compute woutputL by reflecting winputL about m
        pdfOut = specProb;
        if (woutputL.z<DENOM_TOLERANCE) return vec3(0.0);
        float D = microfacetEval(m, roughness);
        float G = smithG2(winputL, woutputL, m, roughness); // Shadow-masking function
        float Fr = fresnelDielectricReflectance(abs(winputL.z), ior);
        vec3 f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(woutputL))*abs(cosTheta(winputL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(winputL, m)), DENOM_TOLERANCE);
        pdfOut *= microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
}

#endif

// ****************************        BSDF common interface        ****************************

vec3 evaluateBsdf( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in int material, in float wavelength_nm, in vec3 rgb, bool fromCamera,
                           inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    evaluateSurface(X, basis, winputL, woutputL, wavelength_nm, rgb); }
#endif
}

vec3 sampleBsdf( in vec3 X, in Basis basis, in vec3 winputL, in int material, in float wavelength_nm, in vec3 rgb, bool fromCamera,
                         inout vec3 woutputL, inout float pdfOut, inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    sampleSurface(X, basis, winputL, wavelength_nm, rgb, woutputL, pdfOut, rnd); }
#endif
}

float pdfBsdf( in vec3 X, in Basis basis, in vec3 winputL, in vec3 woutputL, in int material, in float wavelength_nm, in vec3 rgb, bool fromCamera )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    pdfSurface(X, basis, winputL, woutputL, wavelength_nm, rgb); }
#endif
}

#endif // HAS_GEOMETRY


#ifdef HAS_VOLUME_EMISSION

vec3 VOLUME_EMISSION_EVAL(in vec3 X, in vec3 rgb)
{
    vec3 emission = VOLUME_EMISSION(volumeEmissionColorRGB * volumeEmission, X);
    return emission;
}

#endif // HAS_VOLUME_EMISSION

///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit along specified ray
bool traceRay(in vec3 start, in vec3 dir,
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_INVAL;
    float minMarch = minLengthScale;
    const float HUGE_VAL = 1.0e20;
    float sdf = HUGE_VAL;
#ifdef HAS_SURFACE
    float sdf_surfa = abs(SDF_SURFACE(start));    sdf = min(sdf, sdf_surfa);
#endif
    float InitialSign = sign(sdf);
    float t = InitialSign * sdf; // (always take the first step along the ray direction)
    if (abs(t)>=maxMarchDist) return false;
    for (int n=0; n<__MAX_MARCH_STEPS__; n++)
    {
        vec3 pW = start + t*dir;
        sdf = HUGE_VAL;
#ifdef HAS_SURFACE
        sdf_surfa = abs(SDF_SURFACE(pW));    if (sdf_surfa<minMarch) { material = MAT_SURFA; hit = start + t*dir; return true; } sdf = min(sdf, sdf_surfa);
#endif
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf;
        if (abs(t)>=maxMarchDist) return false;
    }
    return !maxStepsIsMiss;
}

vec3 normal(in vec3 pW, int material)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = 2.0*minLengthScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 Xp = pW+e.xyy; vec3 Xn = pW-e.xyy;
    vec3 Yp = pW+e.yxy; vec3 Yn = pW-e.yxy;
    vec3 Zp = pW+e.yyx; vec3 Zn = pW-e.yyx;
    vec3 N;
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { N = vec3(   SDF_SURFACE(Xp) -    SDF_SURFACE(Xn),    SDF_SURFACE(Yp) -    SDF_SURFACE(Yn),    SDF_SURFACE(Zp) -    SDF_SURFACE(Zn)); return safe_normalize(N); }
#endif
}

#ifdef HAS_NORMALMAP
vec3 perturbNormal(in vec3 X, in Basis basis, int material)
{
#ifdef HAS_SURFACE_NORMALMAP
    if (material==MAT_SURFA)
    {
        vec3 nL_perturbed = SURFACE_NORMAL_MAP(X, basis.nW);
        return localToWorld(normalize(nL_perturbed), basis);
    }
#endif
    return basis.nW;
}
#endif


//////////////////////////////////
// fog transmittance
//////////////////////////////////

// Return the amount of light transmitted (per-channel) along a known "free" segment (i.e. free of geometry)
vec3 transmittanceOverFreeSegment(in vec3 pW, in vec3 rayDir, float segmentLength)
{
    vec3 Tr = vec3(1.0);
#ifdef HAS_FOG
    if (fogEnable)
    {
        vec3 fogAbsorption = vec3(1.0) - fogTint;
        vec3 opticalDepthFog = fogMaxOpticalDepth * (vec3(1.0) - exp(-segmentLength*fogAbsorption/fogMFP));
        Tr *= exp(-opticalDepthFog);
    }
#endif
    return Tr;
}

vec3 transmittanceOverSegment(in vec3 pW, in vec3 rayDir, float segmentLength)
{
    vec3 pW_surface;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_surface, hitMaterial, segmentLength);
    if (hit)
        return vec3(0.0);
    else
#ifdef HAS_FOG
        return transmittanceOverFreeSegment(pW, rayDir, segmentLength);
#else
        return vec3(1.0);
#endif
}

vec3 atmosphericInscatteringRadiance(in vec3 Tr)
{
    if (!fogEnable)
        return vec3(0.0);
    return (1.0 - Tr) * fogEmission;
}



////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
// Sky
//////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dirW)
{
    float rot_phi   = M_PI*envMapPhiRotation/180.0;
    float rot_theta = M_PI*envMapThetaRotation/180.0;
    vec3 sky_pole = vec3(sin(rot_theta), cos(rot_theta), 0.0);
    Basis sky_basis = makeBasis(sky_pole);
    vec3 dirL = worldToLocal(dirW, sky_basis);
    float phi = atan(dirL.y, dirL.x) + M_PI + rot_phi;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = acos(dirL.z);
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 RGB = vec3(1.0);
    if (haveEnvMap)
        RGB = texture(envMap, vec2(u,v)).rgb;
    float t = dot(dirW, sky_pole);
    float tt = envMapTransitionAngle/180.0;
    RGB *= skyPower * mix(skyTintDown, skyTintUp, smoothstep(-tt, tt, t));
    return RGB;
}

vec3 environmentRadiance(in vec3 dir, in vec3 rgb)
{
    vec3 RGB_sky = environmentRadianceRGB(dir);
    return RGB_sky;
}

float skyPowerEstimate()
{
    return 2.0*M_PI * skyPower * 0.5*(maxComponent(skyTintUp) + maxComponent(skyTintDown));
}

vec3 sampleSkyAtSurface(Basis basis, in vec3 rgb, inout vec4 rnd,
                                inout vec3 woutputL, inout vec3 woutputW, inout float pdfDir)
{
    if (skyPower<RADIANCE_EPSILON)
        return vec3(0.0);
    woutputL = sampleHemisphereCosineWeighted(rnd, pdfDir);
    woutputW = localToWorld(woutputL, basis);
    return environmentRadiance(woutputW, rgb);
}


//////////////////////////////////////////////
// Sun
//////////////////////////////////////////////

vec3 sampleSunDir(inout vec4 rnd, inout float pdfDir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    float theta = theta_max * sqrt(rand(rnd));
    float costheta = cos(theta);
    float sintheta = sqrt(max(0.0, 1.0-costheta*costheta));
    float phi = 2.0 * M_PI * rand(rnd);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float x = sintheta * cosphi;
    float y = sintheta * sinphi;
    float z = costheta;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    pdfDir = 1.0/solid_angle;
    return localToWorld(vec3(x, y, z), sunBasis);
}

float pdfSun(in vec3 dir)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return 0.0;
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    return 1.0/solid_angle;
}

vec3 sunRadiance(in vec3 dir, in vec3 rgb)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    if (dot(dir, sunDir) < cos(theta_max)) return vec3(0.0);
    vec3 RGB_sun = sunPower * sunColor;
    return RGB_sun;
}

vec3 sampleSunAtSurface(Basis basis, in vec3 rgb, inout vec4 rnd,
                                inout vec3 woutputL, inout vec3 woutputW, inout float pdfDir)
{
    if (sunPower<RADIANCE_EPSILON)
        return vec3(0.0);
    woutputW = sampleSunDir(rnd, pdfDir);
    woutputL = worldToLocal(woutputW, basis);
    if (woutputL.z < 0.0) return vec3(0.0);
    return sunRadiance(woutputW, rgb);
}


//////////////////////////////////////////////
// Direct lighting routines
//////////////////////////////////////////////

// Estimate direct radiance at the given surface vertex
vec3 directSurfaceLighting(in vec3 pW, Basis basis, in vec3 winputW, in int material,
                                   float wavelength_nm, in vec3 rgb, inout vec4 rnd,
                                   inout float skyPdf, inout float sunPdf, inout float sphPdf)
{
    vec3 winputL = worldToLocal(winputW, basis);
    bool fromCamera = true; // camera path
    vec3 Ldirect = vec3(0.0);
    vec3 woutputL, woutputW;

    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        vec3 Li = sampleSkyAtSurface(basis, rgb, rnd, woutputL, woutputW, skyPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            Li *= transmittanceOverSegment(pW, woutputW, maxLengthScale);
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = max(PDF_EPSILON, pdfBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera));
                vec3 f = evaluateBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera, rnd);
                float misWeight = powerHeuristic(skyPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }
    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        vec3 Li = sampleSunAtSurface(basis, rgb, rnd, woutputL, woutputW, sunPdf);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            Li *= transmittanceOverSegment(pW, woutputW, maxLengthScale);
            if (averageComponent(Li) > RADIANCE_EPSILON)
            {
                // Apply MIS weight with the BSDF pdf for the sampled direction
                float bsdfPdf = max(PDF_EPSILON, pdfBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera));
                vec3 f = evaluateBsdf(pW, basis, winputL, woutputL, material, wavelength_nm, rgb, fromCamera, rnd);
                float misWeight = powerHeuristic(sunPdf, bsdfPdf);
                Ldirect += f * Li/max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
            }
        }
    }
    return min(vec3(radianceClamp), Ldirect);
}


///////////////////////////////////////////////////////////////
// SSS
///////////////////////////////////////////////////////////////

#define HARDMAX_SSS_STEPS 256
#define MIN_SSS_STEPS_BEFORE_RR 4

#ifdef HAS_SURFACE
bool randomwalk_SSS(in vec3 pW, in Basis basis, inout vec4 rnd,
                    in float MFP, in vec3 subsurfaceAlbedo, in vec3 diffuseAlbedoEntry,
                    inout vec3 walk_throughput, inout vec3 pExit)
{
    // Returns true on successful walk to an exit point.
    // Otherwise false indicating termination without finding a valid exit point.
    vec3 pWalk = pW;
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    pWalk -= dPw; // Displace walk into interior to avoid issues with trace numerics
    float pdfDir;
    vec3 dirwalkL = sampleHemisphereCosineWeighted(rnd, pdfDir);
    vec3 dirwalkW = -localToWorld(dirwalkL, basis); // (negative sign to start walk in interior hemisphere)

    // We assume a diffuse Lambertian lobe at the entry point
    // (either pure white, or with the diffuse albedo of the entry point, or somewhere in between)
    vec3 f = mix(vec3(1.0), diffuseAlbedoEntry, subsurfaceDiffuseWeight) / M_PI;
    vec3 fOverPdf = min(vec3(radianceClamp), f/max(PDF_EPSILON, pdfDir));
    vec3 surface_entry_throughput = fOverPdf * abs(dot(dirwalkW, basis.nW));
    walk_throughput = vec3(1.0); //surface_entry_throughput; // update walk throughput due to entry in medium
    for (int n=0; n<HARDMAX_SSS_STEPS; ++n)
    {
        if (n>maxSSSSteps)
            break;
        float walk_step = -log(rand(rnd)) * MFP * lengthScale;
        vec3 pHit;
        int hitMaterial;
        bool hit = traceRay(pWalk, dirwalkW, pHit, hitMaterial, walk_step);
        if (hit)
        {
            // walk left the surface
            pExit = pHit;
            return true;
        }
        // Point remains within the surface, continue walking.
        // First, make a Russian-roulette termination decision (after a minimum number of steps has been taken)
        float termination_prob = 0.0;
        if (n > MIN_SSS_STEPS_BEFORE_RR)
        {
            float continuation_prob = clamp(10.0*maxComponent(walk_throughput), 0.0, 1.0);
            float termination_prob = 1.0 - continuation_prob;
            if (rand(rnd) < termination_prob)
                break;
            walk_throughput /= continuation_prob; // update walk throughput due to RR continuation
        }
        dirwalkW = samplePhaseFunction(dirwalkW, subsurfaceAnisotropy, rnd);
        pWalk += walk_step*dirwalkW;
        walk_throughput *= subsurfaceAlbedo; // update walk throughput due to scattering in medium
    }
    return false;
}

// Estimate radiance at the SSS exit point
vec3 SSS_exit_radiance(in vec3 pW, Basis basis,
                               in vec3 rgb, inout vec4 rnd,
                               in vec3 diffuseAlbedoExit,
                               inout float skyPdf, inout float sunPdf, inout float sphPdf)
{
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    vec3 Ldirect = vec3(0.0);
    vec3 woutputL, woutputW;

    // We assume a diffuse Lambertian lobe at the exit point
    // (either pure white, or with the diffuse albedo of the exit point, or somewhere in between)
    vec3 f = mix(vec3(1.0), diffuseAlbedoExit, subsurfaceDiffuseWeight) / M_PI;

    // Sky
    if (skyPower > RADIANCE_EPSILON)
    {
        vec3 Li = sampleSkyAtSurface(basis, rgb, rnd, woutputL, woutputW, skyPdf);
        Li *= transmittanceOverSegment(pW, woutputW, maxLengthScale);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            // Apply MIS weight with the BSDF pdf for the sampled direction
            float bsdfPdf = pdfHemisphereCosineWeighted(woutputL);
            float misWeight = powerHeuristic(skyPdf, bsdfPdf);
            Ldirect += f * Li/max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
        }
    }
    // Sun
    if (sunPower > RADIANCE_EPSILON)
    {
        vec3 Li = sampleSunAtSurface(basis, rgb, rnd, woutputL, woutputW, sunPdf);
        Li *= transmittanceOverSegment(pW, woutputW, maxLengthScale);
        if (averageComponent(Li) > RADIANCE_EPSILON)
        {
            // Apply MIS weight with the BSDF pdf for the sampled direction
            float bsdfPdf = pdfHemisphereCosineWeighted(woutputL);
            float misWeight = powerHeuristic(sunPdf, bsdfPdf);
            Ldirect += f * Li/max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, basis.nW)) * misWeight;
        }
    }
    return min(vec3(radianceClamp), Ldirect);
}

#endif


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
    primaryDir = safe_normalize(camDir + s);
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
    primaryDir = safe_normalize(focalPlaneHit - lensPos);
}


vec3 cameraPath(in vec3 primaryStart, in vec3 primaryDir, 
                        float wavelength_nm, in vec3 rgb, inout vec4 rnd)
{
    // Perform pathtrace starting from the camera lens to estimate the primary ray radiance, L
    vec3 L = vec3(0.0);
    float misWeightSky = 1.0; // For MIS book-keeping
    float misWeightSun = 1.0; // For MIS book-keeping
    vec3 pW = primaryStart;
    vec3 rayDir = primaryDir; // (opposite to light beam direction)

    vec3 throughput = vec3(1.0);
    for (int vertex=0; vertex<=__MAX_BOUNCES__; ++vertex)
    {
        // Bail out now if the path continuation throughput is below threshold
        if (maxComponent(throughput) < THROUGHPUT_EPSILON) break;

        // Raycast along current propagation direction rayDir, from current vertex pW to pW_next
        vec3 pW_next;
        int hitMaterial;
        bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
        float rayLength = maxLengthScale;

        vec3 Tr;
        if (hit)
        {
            rayLength = length(pW_next - pW);
            Tr = transmittanceOverFreeSegment(pW, rayDir, rayLength);
        }
        else
        {
            Tr = transmittanceOverFreeSegment(pW, rayDir, maxLengthScale);
#ifdef HAS_FOG
            if (averageComponent(Tr) > RADIANCE_EPSILON)
#endif
            {
                // This ray missed all geometry; add contribution from distant lights and terminate path
                if (!(vertex==0 && !envMapVisible)      && misWeightSky>0.0) L += throughput * Tr * misWeightSky * environmentRadiance(rayDir, rgb);
                if (!(vertex==0 && !sunVisibleDirectly) && misWeightSun>0.0) L += throughput * Tr * misWeightSun * sunRadiance(rayDir, rgb);
            }
#ifdef HAS_FOG
            // Add term for in-scattering in the homogeneous atmosphere (fog) over the segment to the light hit
            if (vertex<2) // (restrict to first two segments only, for efficiency)
                L += throughput * atmosphericInscatteringRadiance(Tr);
#endif
            // Ray escapes to infinity
            break;
        }

#ifdef HAS_FOG
        // Add term for in-scattering in the homogeneous atmosphere (fog) over the segment to the light hit
        if (vertex<2) // (restrict to first two segments only, for efficiency)
            L += throughput * atmosphericInscatteringRadiance(Tr);

        // Attenuate throughput to next hit due to atmospheric attenuation (Beer's law)
        throughput *= Tr;
#endif

        // This ray hit some geometry, so deal with the surface interaction.
        // First, compute the normal and thus the local vertex basis:
        pW = pW_next;
        vec3 nW = normal(pW, hitMaterial);
        vec3 ngW = nW; // geometric normal, needed for robust ray continuation
        Basis basis = makeBasis(nW);
#ifdef HAS_NORMALMAP
        nW = perturbNormal(pW, basis, hitMaterial);
        // If the incident ray lies below the hemisphere of the perturbed shading normal,
        // which can occur due to normal mapping, apply the "Flipping hack" to prevent artifacts
        // (see Schüßler, "Microfacet-based Normal Mapping for Robust Monte Carlo Path Tracing")
        if ((dot(nW, rayDir) > 0.0) && (hitMaterial != MAT_DIELE))
            nW = 2.0*ngW*dot(ngW, nW) - nW;
        basis = makeBasis(nW);
#endif

        // Make a binary choice whether to scatter at the surface, or do a subsurface random walk:
        bool do_subsurface_walk = false;
        float prob_sss = 0.0;
        vec3 subsurfaceAlbedo;
#ifdef HAS_SURFACE
        if (subsurfaceMFP > 0.0)         // a) the surface must have non-zero subsurface MFP, and
        {
            subsurfaceAlbedo = SUBSURFACE_ALBEDO_EVAL(pW, basis.nW, rgb);
            float diffuse_weight    = (1.0 - subsurface)*averageComponent(SURFACE_DIFFUSE_REFL_EVAL(pW, basis.nW, -rayDir, rgb));
            float    spec_weight    = averageComponent(SURFACE_SPEC_REFL_EVAL(pW, basis.nW, -rayDir, rgb));
            float subsurface_weight = subsurface * averageComponent(subsurfaceAlbedo) * 3.0; // weight more highly due to higher variance
            float total_weight = diffuse_weight + spec_weight + subsurface_weight;
            prob_sss = subsurface_weight / (total_weight + DENOM_TOLERANCE);
            prob_sss = clamp(prob_sss, 0.0, 1.0-PDF_EPSILON);
            do_subsurface_walk = (rand(rnd) < prob_sss);
        }
#endif

        // Regular surface bounce case
        if (!do_subsurface_walk)
        {
            // Sample BSDF for the next bounce direction
            vec3 winputW = -rayDir; // winputW, points *towards* the incident direction
            vec3 winputL = worldToLocal(winputW, basis);
            vec3 woutputL; // woutputL, points *towards* the outgoing direction
            bool fromCamera = true; // camera path
            float bsdfPdf;
            vec3 f = sampleBsdf(pW, basis, winputL, hitMaterial, wavelength_nm, rgb, fromCamera, woutputL, bsdfPdf, rnd);
            vec3 woutputW = localToWorld(woutputL, basis);

            // Detect dielectric transmission
#ifdef HAS_VOLUME_EMISSION
            // Add volumetric emission at the surface point, if present (treating it as an isotropic radiance field)
            L += throughput * VOLUME_EMISSION_EVAL(pW, rgb);
#endif

            // Update ray direction to the BSDF-sampled direction
            rayDir = woutputW;

            // Prepare for tracing the direct lighting and continuation rays
            pW += ngW * sign(dot(rayDir, ngW)) * 3.0*minLengthScale; // perturb vertex into half-space of scattered ray

            // Add direct lighting term at current surface vertex
            float skyPdf = 0.0;
            float sunPdf = 0.0;
            float sphPdf = 0.0;
            L += throughput * directSurfaceLighting(pW, basis, winputW, hitMaterial, wavelength_nm, rgb, rnd,
                                                    skyPdf, sunPdf, sphPdf);

            // Update path continuation throughput
            vec3 fOverPdf = min(vec3(radianceClamp), f/max(PDF_EPSILON, bsdfPdf));
            vec3 surface_throughput = fOverPdf * abs(dot(woutputW, nW)) / max(PDF_EPSILON, 1.0 - prob_sss);
            throughput *= surface_throughput;

            misWeightSky = powerHeuristic(bsdfPdf, skyPdf); // compute sky MIS weight for bounce ray
            misWeightSun = powerHeuristic(bsdfPdf, sunPdf); // compute sun MIS weight for bounce ray
        }

#ifdef HAS_SURFACE
        // Do subsurface random walk to a new vertex
        else
        {
            vec3 diffuseAlbedoEntry = SURFACE_DIFFUSE_REFL_EVAL(pW, basis.nW, -rayDir, rgb);
            vec3 walk_throughput;
            vec3 pExit;
            basis.nW = ngW; // (use basis with geometric normal for SSS entry, to avoid surface acne)
            bool success = randomwalk_SSS(pW, basis, rnd, subsurfaceMFP, subsurfaceAlbedo, diffuseAlbedoEntry, walk_throughput, pExit);
            if (!success)
                break;

            // Update path vertex to exit point
            pW = pExit;

            // Compute updated normal and basis at exit point
            nW = normal(pExit, MAT_SURFA);
            Basis basis_exit = makeBasis(nW);

            // Sample direction of exit ray (assuming a diffuse Lambertian lobe with diffuse albedo of exit point)
            float bsdfPdf;
            vec3 dirExitL = sampleHemisphereCosineWeighted(rnd, bsdfPdf);
            vec3 woutputW = localToWorld(dirExitL, basis_exit);

            // Update ray direction to the SSS exit direction
            rayDir = woutputW;

            // Update throughput due to random walk
            throughput *= walk_throughput;

            // Add direct lighting term at exit vertex (assumed to be a diffuse lobe)
            float skyPdf = 0.0;
            float sunPdf = 0.0;
            float sphPdf = 0.0;
            vec3 diffuseAlbedoExit = SURFACE_DIFFUSE_REFL_EVAL(pW, nW, -rayDir, rgb);
            L += throughput * SSS_exit_radiance(pExit, basis_exit, rgb, rnd, diffuseAlbedoExit,
                                                skyPdf, sunPdf, sphPdf);
            misWeightSky = powerHeuristic(bsdfPdf, skyPdf); // compute sky MIS weight for bounce ray
            misWeightSun = powerHeuristic(bsdfPdf, sunPdf); // compute sun MIS weight for bounce ray

            // Update path continuation throughput
            vec3 f = mix(vec3(1.0), diffuseAlbedoExit, subsurfaceDiffuseWeight) / M_PI;
            vec3 fOverPdf = min(vec3(radianceClamp), f/max(PDF_EPSILON, bsdfPdf));
            vec3 surface_exit_throughput = fOverPdf * abs(dot(woutputW, nW));
            throughput *= surface_exit_throughput / max(PDF_EPSILON, prob_sss);

            // Prepare for tracing the continuation ray
            pW += nW * sign(dot(rayDir, nW)) * 3.0*minLengthScale; // perturb vertex into half-space of scattered ray
        }
#endif
    }
    return L;
}

float sample_jitter(float xi)
{
    // sample from triangle filter, returning jitter in [-1, 1]
    return xi < 0.5 ? sqrt(2.0*xi) - 1.0 : 1.0 - sqrt(2.0 - 2.0*xi);
}

void pathtrace(vec2 pixel, vec4 rnd) // the current pixel
{
    float wavelength_nm = 390.0 + (750.0 - 390.0)*0.5; // non-dispersive rendering uses mid-range wavelength
    vec3 rgb = vec3(0.0);

    // Setup sun basis
    sunBasis = makeBasis(sunDir);

    // Sample radiance of primary ray
    vec3 L = vec3(0.0);
    for (int n=0; n<__MAX_SAMPLES_PER_FRAME__; ++n)
    {
        // Apply FIS to obtain pixel jitter about center in pixel units
        float jx = 0.5 * filterRadius * sample_jitter(rand(rnd));
        float jy = 0.5 * filterRadius * sample_jitter(rand(rnd));
        vec2 pixelj = pixel + vec2(jx, jy);

        // Compute world ray direction for this fragment
        vec3 primaryStart, primaryDir;
#ifdef HAS_CUSTOM_CAMERA
        CONSTRUCT_PRIMARY_RAY(pixelj, rnd, primaryStart, primaryDir);
#else
        constructPrimaryRay(pixelj, rnd, primaryStart, primaryDir);
#endif

        // Perform pathtrace to estimate the primary ray radiance, L
        L += cameraPath(primaryStart, primaryDir, wavelength_nm, rgb, rnd);
    }
    L /= float(__MAX_SAMPLES_PER_FRAME__);

    // Write updated radiance and sample count
    vec4 oldL = texture(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;

    // Compute tristimulus contribution from estimated radiance
    vec3 colorXYZ = rgbToXyz(L);
    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;
    gbuf_rad = vec4(newL, newN);
    gbuf_rng = rnd;
}

void main()
{
    vec4 rnd = texture(RngData, vTexCoord);
    INIT();
    pathtrace(gl_FragCoord.xy, rnd);
}
`,

'simplepathtracer-vertex-shader': `#version 300 es
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
uniform float contrast;
uniform float saturation;
uniform float hueShift;

out vec4 g_outputColor;

float toneMap(float L)
{
  return L / (1.0 + L);
}

vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));
    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{
    vec3 L = texture(Radiance, vTexCoord).rgb;
    float X = L.x;
    float Y = L.y;
    float Z = L.z;
    
    // convert XYZ tristimulus to linear RGB color space
    vec3 RGB;
    RGB.r =  3.2406*X - 1.5372*Y - 0.4986*Z;
    RGB.g = -0.9689*X + 1.8758*Y + 0.0415*Z;
    RGB.b =  0.0557*X - 0.2040*Y + 1.0570*Z;
    
    // deal with out-of-gamut RGB.
    float delta = -min(0.0, min(min(RGB.r, RGB.g), RGB.b));
    RGB.r += delta;
    RGB.g += delta;
    RGB.b += delta;

    // apply gamma correction to convert linear RGB to sRGB
    RGB = pow(RGB, vec3(invGamma));

    // apply tonemapping
    RGB *= pow(2.0, exposure);
    float R = RGB.r;
    float G = RGB.g;
    float B = RGB.b;
    R = toneMap(R);
    G = toneMap(G);
    B = toneMap(B);



    // apply saturation
    float mean = (R + G + B)/3.0;
    float dR = R - mean;
    float dG = G - mean;
    float dB = B - mean;
    R = mean + sign(dR)*pow(abs(dR), 1.0/saturation);
    G = mean + sign(dG)*pow(abs(dG), 1.0/saturation);
    B = mean + sign(dB)*pow(abs(dB), 1.0/saturation);

    // apply contrast
    dR = R - 0.5;
    dG = G - 0.5;
    dB = B - 0.5;
    R = 0.5 + sign(dR)*pow(abs(dR), 1.0/contrast);
    G = 0.5 + sign(dG)*pow(abs(dG), 1.0/contrast);
    B = 0.5 + sign(dB)*pow(abs(dB), 1.0/contrast);

    vec3 C =vec3(R,G,B);

    // apply hue-shift
    if (hueShift > 0.0)
    {
        vec3 hsv = rgb2hsv(C);
        hsv.r = mod(hsv.r + hueShift, 1.0);
        C = hsv2rgb(hsv);
    }

    g_outputColor = vec4(C, 1.0);
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