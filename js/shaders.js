var Shaders = {

'pathtracer-fragment-shader': `

#extension GL_EXT_draw_buffers : require
precision highp float;

uniform sampler2D Radiance;         // 0
uniform sampler2D RngData;          // 1
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
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
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

IOR_FUNC

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

float pdfHemisphere(in vec3 wiL)
{
    // pdf for cosine-weighted sampling of hemisphere
    return wiL.z / M_PI;
}

float powerHeuristic(const float a, const float b)
{
    float t = a;
    return t / (t + b);
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

float evaluateDielectric( in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    float Fr = DIELECTRIC_FRESNEL(X, nW) * fresnelDielectricReflectance(woL.z, ior, 1.0);
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
    float roughness = dieleRoughness * DIELECTRIC_ROUGHNESS(X, nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f;
    if (reflected)
    {
        f = Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
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

float pdfDielectric( in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    float Fr = DIELECTRIC_FRESNEL(X, nW) * fresnelDielectricReflectance(woL.z, ior, 1.0);
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
    float roughness = dieleRoughness * DIELECTRIC_ROUGHNESS(X, nW);
    pdf *= microfacetPDF(h, roughness);
    return abs(pdf * dwh_dwo);
}

float sampleDielectric( in vec3 X, in vec3 nW, in vec3 woL, in float wavelength_nm, in vec3 XYZ,
                        inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_DIELE(wavelength_nm);
    float Fr = DIELECTRIC_FRESNEL(X, nW) * fresnelDielectricReflectance(woL.z, ior, 1.0);
    float roughness = dieleRoughness * DIELECTRIC_ROUGHNESS(X, nW);
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

float evaluateMetal( in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = METAL_FRESNEL(X, nW) * fresnelMetalReflectance(woL.z, ior, k);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float roughness = metalRoughness * METAL_ROUGHNESS(X, nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f = Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    return f;
}

float pdfMetal( in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 XYZ )
{
    float ior = IOR_DIELE(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float roughness = metalRoughness * METAL_ROUGHNESS(X, nW);
    float pdf = microfacetPDF(h, roughness) * dwh_dwo;
    return pdf;
}

float sampleMetal( in vec3 X, in vec3 nW, in vec3 woL, in float wavelength_nm, in vec3 XYZ,
                   inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = METAL_FRESNEL(X, nW) * fresnelMetalReflectance(woL.z, ior, k);
    float roughness = metalRoughness * METAL_ROUGHNESS(X, nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
    if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
    float D = microfacetEval(m, roughness);
    float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
    float f = Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    float dwh_dwo; // Jacobian of the half-direction mapping
    dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
    pdfOut = microfacetPDF(m, roughness) * dwh_dwo;
    return f;
}


// ****************************        Surface        ****************************

vec3 SURFACE_DIFFUSE_REFL_XYZ(in vec3 X, in vec3 nW)
{
    vec3 reflRGB = xyzToRgb(surfaceDiffuseAlbedoXYZ) * SURFACE_DIFFUSE_REFLECTANCE(X, nW);
    vec3 reflXYZ = rgbToXyz(reflRGB);
    return xyz_to_spectrum(reflXYZ);
}

vec3 SURFACE_SPEC_REFL_XYZ(in vec3 X, in vec3 nW)
{
    vec3 reflRGB = xyzToRgb(surfaceSpecAlbedoXYZ) * SURFACE_SPECULAR_REFLECTANCE(X, nW);
    vec3 reflXYZ = rgbToXyz(reflRGB);
    return xyz_to_spectrum(reflXYZ);
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

float evaluateSurface(in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    float diffuseAlbedo = clamp(dot(XYZ, SURFACE_DIFFUSE_REFL_XYZ(X, nW)), 0.0, 1.0);
    float    specAlbedo = clamp(dot(XYZ, SURFACE_SPEC_REFL_XYZ(X, nW)),    0.0, 1.0);
    float ior = surfaceIor;
    float roughness = surfaceRoughness * SURFACE_ROUGHNESS(X, nW);
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float specWeight = specAlbedo * Fr;
    float diffWeight = diffuseAlbedo * (1.0-specWeight);
    float f = specWeight * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    f += diffWeight/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    float diffuseAlbedo = clamp(dot(XYZ, SURFACE_DIFFUSE_REFL_XYZ(X, nW)), 0.0, 1.0);
    float    specAlbedo = clamp(dot(XYZ, SURFACE_SPEC_REFL_XYZ(X, nW)),    0.0, 1.0);
    float ior = surfaceIor;
    float roughness = surfaceRoughness * SURFACE_ROUGHNESS(X, nW);
    float diffusePdf = pdfHemisphere(wiL);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    float specWeight = specAlbedo * Fr;
    float diffWeight = diffuseAlbedo * (1.0-specWeight);
    float weightSum = max(specWeight + diffWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

float sampleSurface(in vec3 X, in vec3 nW, in vec3 woL, in vec3 XYZ,
                    inout vec3 wiL, inout float pdfOut, inout vec4 rnd)
{
    float diffuseAlbedo = clamp(dot(XYZ, SURFACE_DIFFUSE_REFL_XYZ(X, nW)), 0.0, 1.0);
    float    specAlbedo = clamp(dot(XYZ, SURFACE_SPEC_REFL_XYZ(X, nW)),    0.0, 1.0);
    float ior = surfaceIor;
    float roughness = surfaceRoughness * SURFACE_ROUGHNESS(X, nW);
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    float specWeight = specAlbedo * Fr;
    float diffWeight = diffuseAlbedo * (1.0-specWeight);
    float weightSum = max(specWeight + diffWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        wiL = sampleHemisphere(rnd, pdfOut);
        pdfOut *= (1.0-specProb);
        return diffWeight/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
        if (wiL.z<DENOM_TOLERANCE) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float f = specWeight * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = specProb * microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
}

// ****************************        BSDF common interface        ****************************

float evaluateBsdf( in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ, 
                    inout vec4 rnd )
{
    if      (material==MAT_DIELE) { return evaluateDielectric(X, nW, woL, wiL, wavelength_nm, XYZ); }
    else if (material==MAT_METAL) { return      evaluateMetal(X, nW, woL, wiL, wavelength_nm, XYZ); }
    else                          { return    evaluateSurface(X, nW, woL, wiL,                XYZ); }
}

float sampleBsdf( in vec3 X, in vec3 nW, in vec3 woL, in int material, in float wavelength_nm, in vec3 XYZ,
                  inout vec3 wiL, inout float pdfOut, inout vec4 rnd ) 
{
    if      (material==MAT_DIELE) { return sampleDielectric(X, nW, woL, wavelength_nm, XYZ, wiL, pdfOut, rnd); }
    else if (material==MAT_METAL) { return      sampleMetal(X, nW, woL, wavelength_nm, XYZ, wiL, pdfOut, rnd); }
    else                          { return    sampleSurface(X, nW, woL,                XYZ, wiL, pdfOut, rnd); }
}

float pdfBsdf( in vec3 X, in vec3 nW, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ )
{
    if      (material==MAT_DIELE) { return pdfDielectric(X, nW, woL, wiL, wavelength_nm, XYZ); }
    else if (material==MAT_METAL) { return      pdfMetal(X, nW, woL, wiL, wavelength_nm, XYZ); }
    else                          { return    pdfSurface(X, nW, woL, wiL,                XYZ); }
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

float environmentRadiance(in vec3 dir, in vec3 XYZ)
{
    vec3 XYZ_sky = environmentRadianceXYZ(dir);

    // convert to radiance at the given wavelength
    vec3 XYZ_spec = xyz_to_spectrum(XYZ_sky);
    return dot(XYZ, XYZ_spec);
}

float directLighting(in vec3 pW, Basis basis, in vec3 woW, in int material, 
                     float wavelength_nm, in vec3 XYZ, inout vec4 rnd)
{
    float lightPdf;
    float Li;
    vec3 wiW; // direction of sampled direct light (*towards* the light)
    {
        // Env-map sampling
        float hemispherePdf;
        vec3 wiL = sampleHemisphere(rnd, hemispherePdf);
        lightPdf = hemispherePdf;
        wiW = localToWorld(wiL, basis);
        Li = environmentRadiance(wiW, XYZ);

        bool occluded = Occluded(pW, wiW); 
        if (occluded) Li *= abs(1.0 - shadowStrength);
    }

    // Apply MIS weight with the BSDF pdf for the sampled direction
    vec3 woL = worldToLocal(woW, basis);
    vec3 wiL = worldToLocal(wiW, basis);
    float bsdfPdf = pdfBsdf(pW, basis.nW, woL, wiL, material, wavelength_nm, XYZ);
    if ( bsdfPdf<PDF_EPSILON ) return 0.0;

    float f = evaluateBsdf(pW, basis.nW, woL, wiL, material, wavelength_nm, XYZ, rnd);
    float misWeight = powerHeuristic(lightPdf, bsdfPdf);
    float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li * abs(dot(wiW, basis.nW)) * misWeight;
}

////////////////////////////////////////////////////////////////////////////////
// Pathtracing integrator
////////////////////////////////////////////////////////////////////////////////

vec3 constuctPrimaryDir(in vec2 pixel)
{
    // Compute world ray direction for given (possibly jittered) fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
    float fh = camNear*tan(0.5*radians(camFovy)) / camZoom; // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    return normalize(camNear*camDir + s);
}

void pathtrace(vec2 pixel, vec4 rnd) // the current pixel
{
    if (rand(rnd) < skipProbability)
    {
        vec4 oldL = texture2D(Radiance, vTexCoord);
        float oldN = oldL.w;
        float newN = oldN;
        vec3 newL = oldL.rgb;

        gl_FragData[0] = vec4(newL, newN);
        gl_FragData[1] = rnd;
        return;
    } 

    // Sample photon wavelength via the inverse CDF of the emission spectrum
    // (here w is spectral offset, i.e. wavelength = 360.0 + (750.0 - 360.0)*w)
    // (linear interpolation into the inverse CDF texture and XYZ table should ensure smooth sampling over the range)
    float w = texture2D(ICDF, vec2(rand(rnd), 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;

    // Convert wavelength to XYZ tristimulus
    vec3 XYZ = texture2D(WavelengthToXYZ, vec2(w, 0.5)).rgb;

    // Jitter over pixel
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));

    // Compute world ray direction for this fragment
    vec3 primaryDir = constuctPrimaryDir(pixel);

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int rayMaterial = MAT_VACUU;
    int hitMaterial;
    bool hit = traceRay(camPos, primaryDir, pW, hitMaterial, maxScale);
    
    vec3 colorXYZ; 
    if ( !hit )
    {
        if (envMapVisible) colorXYZ = environmentRadianceXYZ(primaryDir);
        else               colorXYZ = vec3(0.0);
    }
    else
    {   
        float L = 0.0;    
        float throughput = 1.0;
        for (int bounce=0; bounce<MAX_BOUNCES; ++bounce)
        {
            // Compute normal at current surface vertex
            vec3 nW = normal(pW, hitMaterial);
            Basis basis = makeBasis(nW);

            // Add direct lighting term at current surface vertex
            L += throughput * directLighting(pW, basis, woW, hitMaterial, wavelength_nm, XYZ, rnd);

            // Sample BSDF for the next bounce direction
            vec3 woL = worldToLocal(woW, basis);
            vec3 wiL;
            float bsdfPdf;
            float f = sampleBsdf(pW, nW, woL, hitMaterial, wavelength_nm, XYZ, wiL, bsdfPdf, rnd);
            vec3 wiW = localToWorld(wiL, basis);

            // material in which ray propagates changes (only) on transmission
            if (dot(woW, nW) * dot(wiW, nW) < 0.0) 
            {
                rayMaterial = hitMaterial;
            }

            // Update path throughput
            float fOverPdf = min(radianceClamp, f/max(PDF_EPSILON, bsdfPdf));
            throughput *= fOverPdf * abs(dot(wiW, nW));
            if (throughput < THROUGHPUT_EPSILON) break;

            // Trace bounce ray
            float displacement = 3.0*minScale;
            float wiWnW = dot(wiW, nW);
            pW += nW * sign(wiWnW) * displacement; // perturb vertex into half-space of scattered ray
            vec3 pW_next;
            bool hit = traceRay(pW, wiW, pW_next, hitMaterial, maxScale);

            // If ray missed, add environment light term and terminate path
            if (!hit)
            {
                float lightPdf = pdfHemisphere(wiL);
                float misWeight = powerHeuristic(bsdfPdf, lightPdf);
                float Li = environmentRadiance(wiW, XYZ);
                L += throughput * Li * misWeight;
                break;
            }

            // Support non-physical illumination of occluded areas by the env map
            else if (shadowStrength < 1.0)
            {
                float lightPdf = pdfHemisphere(wiL);
                float misWeight = powerHeuristic(bsdfPdf, lightPdf);
                float Li = environmentRadiance(wiW, XYZ);
                L += throughput * Li * misWeight * abs(1.0 - shadowStrength);
            }

            // If the bounce ray lies inside a dielectric, apply Beer's law for absorption
            if (rayMaterial==MAT_DIELE)
            {
                float absorptionLength = length(pW_next - pW);
                vec3 RGB = clamp(xyzToRgb(XYZ), 0.0, 1.0);
                throughput *= exp(-absorptionLength*dot(dieleAbsorptionRGB, RGB));
            }

            // Update vertex
            vec3 rayDir = normalize(pW_next - pW);
            woW = -rayDir;
            pW = pW_next;
        }

        colorXYZ = XYZ * L;
    }

    // Write updated radiance and sample count
    vec4 oldL = texture2D(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + colorXYZ) / newN;

    gl_FragData[0] = vec4(newL, newN);
    gl_FragData[1] = rnd;
}


void ENTRY_PATHTRACE_ALL()
{
    vec4 rnd = texture2D(RngData, vTexCoord);
    pathtrace(gl_FragCoord.xy, rnd);
}


////////////////////////////////////////////////////////////////////////////////
// Normals integrator
////////////////////////////////////////////////////////////////////////////////

void ENTRY_NORMALS()
{
    vec4 rnd = texture2D(RngData, vTexCoord);
    vec2 pixel = gl_FragCoord.xy;
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));
    vec3 primaryDir = constuctPrimaryDir(pixel);

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int rayMaterial = MAT_VACUU;
    int hitMaterial;
    bool hit = traceRay(camPos, primaryDir, pW, hitMaterial, maxScale);
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


////////////////////////////////////////////////////////////////////////////////
// Ambient occlusion integrator
////////////////////////////////////////////////////////////////////////////////

void ENTRY_AO()
{
    vec4 rnd = texture2D(RngData, vTexCoord);
    if (rand(rnd) < skipProbability)
    {
        vec4 oldL = texture2D(Radiance, vTexCoord);
        float oldN = oldL.w;
        float newN = oldN;
        vec3 newL = oldL.rgb;

        gl_FragData[0] = vec4(newL, newN);
        gl_FragData[1] = rnd;
        return;
    } 

    // Jitter over pixel
    vec2 pixel = gl_FragCoord.xy;
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));
    vec3 primaryDir = constuctPrimaryDir(pixel);
    
    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int rayMaterial = MAT_VACUU;
    int hitMaterial;
    bool hit = traceRay(camPos, primaryDir, pW, hitMaterial, maxScale);

    vec3 colorXYZ;
    if (hit)
    {
        // Compute normal at hit point
        vec3 nW = normal(pW, hitMaterial);
        Basis basis = makeBasis(nW);

        // Construct a uniformly sampled AO ray
        float hemispherePdf;
        vec3 wiL = sampleHemisphere(rnd, hemispherePdf);
        vec3 wiW = localToWorld(wiL, basis);

        // Compute diffuse albedo at sampled wavelength
        float w = texture2D(ICDF, vec2(rand(rnd), 0.5)).r;
        vec3 XYZ = texture2D(WavelengthToXYZ, vec2(w, 0.5)).rgb;
        float diffuseAlbedo = clamp(dot(XYZ, SURFACE_DIFFUSE_REFL_XYZ(pW, nW)), 0.0, 1.0);

        // Set incident radiance to according to whether the AO ray hit anything or missed.
        float L;
        if (!Occluded(pW, wiW)) L = diffuseAlbedo;
        else                    L = diffuseAlbedo * abs(1.0 - shadowStrength);
        colorXYZ = XYZ * L;
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


////////////////////////////////////////////////////////////////////////////////
// BRDF integrator
////////////////////////////////////////////////////////////////////////////////

void ENTRY_BRDF()
{
    vec4 rnd = texture2D(RngData, vTexCoord);
    if (rand(rnd) < skipProbability)
    {
        vec4 oldL = texture2D(Radiance, vTexCoord);
        float oldN = oldL.w;
        float newN = oldN;
        vec3 newL = oldL.rgb;

        gl_FragData[0] = vec4(newL, newN);
        gl_FragData[1] = rnd;
        return;
    } 

    // Jitter over pixel
    vec2 pixel = gl_FragCoord.xy;
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));
    vec3 primaryDir = constuctPrimaryDir(pixel);

    float w = texture2D(ICDF, vec2(rand(rnd), 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;
    vec3 XYZ = texture2D(WavelengthToXYZ, vec2(w, 0.5)).rgb;
    
    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int rayMaterial = MAT_VACUU;
    int hitMaterial;
    bool hit = traceRay(camPos, primaryDir, pW, hitMaterial, maxScale);

    vec3 colorXYZ;
    if (hit)
    {
        // Compute normal at hit point
        vec3 nW = normal(pW, hitMaterial);
        Basis basis = makeBasis(nW);

        vec3 woL = worldToLocal(woW, basis);
        vec3 wiL;
        float bsdfPdf;
        float f = sampleBsdf(pW, nW, woL, hitMaterial, wavelength_nm, XYZ, wiL, bsdfPdf, rnd);
        colorXYZ = XYZ * f;
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


////////////////////////////////////////////////////////////////////////////////
// Surface diffuse integrator
////////////////////////////////////////////////////////////////////////////////

void ENTRY_SURFACE_DIFFUSE()
{
    vec4 rnd = texture2D(RngData, vTexCoord);
    if (rand(rnd) < skipProbability)
    {
        vec4 oldL = texture2D(Radiance, vTexCoord);
        float oldN = oldL.w;
        float newN = oldN;
        vec3 newL = oldL.rgb;

        gl_FragData[0] = vec4(newL, newN);
        gl_FragData[1] = rnd;
        return;
    } 

    // Jitter over pixel
    vec2 pixel = gl_FragCoord.xy;
    if (jitter) pixel += (-0.5 + vec2(rand(rnd), rand(rnd)));
    vec3 primaryDir = constuctPrimaryDir(pixel);

    float w = texture2D(ICDF, vec2(rand(rnd), 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;
    vec3 XYZ = texture2D(WavelengthToXYZ, vec2(w, 0.5)).rgb;
    
    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int rayMaterial = MAT_VACUU;
    int hitMaterial;
    bool hit = traceRay(camPos, primaryDir, pW, hitMaterial, maxScale);

    vec3 colorXYZ;
    if (hit)
    {
        vec3 nW = normal(pW, hitMaterial);
        colorXYZ = SURFACE_DIFFUSE_REFL_XYZ(pW, nW);
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
`,

'pathtracer-vertex-shader': `

#extension GL_EXT_draw_buffers : require
precision highp float;

attribute vec3 Position;
attribute vec2 TexCoord;

varying vec2 vTexCoord;

void main() 
{
	gl_Position = vec4(Position, 1.0);
	vTexCoord = TexCoord;
}
`,

'tonemapper-fragment-shader': `

#extension GL_EXT_draw_buffers : require
precision highp float;

uniform sampler2D Radiance;
varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float whitepoint;


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
	vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
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

	gl_FragColor = vec4(S, 0.0);
}
`,

'tonemapper-vertex-shader': `

#extension GL_EXT_draw_buffers : require
precision highp float;

attribute vec3 Position;
attribute vec2 TexCoord;
varying vec2 vTexCoord;

void main() 
{
	gl_Position = vec4(Position, 1.0);
	vTexCoord = TexCoord;
}
`,

}