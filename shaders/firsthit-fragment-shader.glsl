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

uniform vec3 surfaceDiffuseAlbedoRGB;
uniform vec3 surfaceSpecAlbedoRGB;
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
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, woW);
    return reflRGB;
}

vec3 SURFACE_SPEC_REFL_RGB(in vec3 X, in vec3 nW, in vec3 woW)
{
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(surfaceSpecAlbedoRGB, X, nW, woW);
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
    pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

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
