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

#ifdef HAS_VOLUME_EMISSION
vec3 VOLUME_EMISSION_EVAL(in vec3 X)
{
    vec3 emission = VOLUME_EMISSION(volumeEmissionColorRGB * volumeEmission, X);
    return emission;
}
#endif // HAS_VOLUME_EMISSION

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
    vec3 Tr_pW;
    float invDistancePdf;
    if (T01 > 1.0e-3)
    {
        float t_scatter = t0 - log(1.0 - rand(rnd)*(1.0 - T01)) / extinction_norm; // sampled scatter distance
        vec3 pW_scatter = pW + t_scatter*rayDir;                                   // sampled scatter point
        vec3 opticalDepth_scatter = (t_scatter - t0) * extinction;                 // optical depth from pW -> pW_scatter
        Tr_pW = exp(-opticalDepth_scatter);                                        // transmittance from pW -> pW_scatter
        invDistancePdf = (1.0 - T01) * exp(averageComponent(opticalDepth_scatter)) / extinction_norm; // PDF of scatter distance
    }
    else
    {
        
    }

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
#ifdef HAS_ATMOSPHERE
                vec3 TrToLight = transmittanceOverSegment(pW_hit+dPw, woutputW, maxLengthScale);
#else
                vec3 TrToLight = vec3(1.0);
#endif
                Ldirect += f * TrToLight * Li / max(PDF_EPSILON, skyPdf) * abs(dot(woutputW, nW));
            }
            if (sunPower > RADIANCE_EPSILON)
            {
                // Sun
                float sunPdf;
                vec3 woutputW;
                vec3 Li = sampleSunAtSurface(basis, rnd, woutputW, sunPdf);
#ifdef HAS_ATMOSPHERE
                vec3 TrToLight = transmittanceOverSegment(pW_hit+dPw, woutputW, maxLengthScale);
#else
                vec3 TrToLight = vec3(1.0);
#endif
                if (averageComponent(Li) > RADIANCE_EPSILON)
                {
                    Li *= Transmittance(pW_hit+dPw, woutputW);
                    Ldirect += f * TrToLight * Li / max(PDF_EPSILON, sunPdf) * abs(dot(woutputW, nW));
                }
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
