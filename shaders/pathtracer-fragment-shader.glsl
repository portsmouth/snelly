precision highp float;

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
uniform float envMapPhiRotation;
uniform float envMapThetaRotation;
uniform float envMapTransitionAngle;
uniform float skyPower;
uniform vec3 skyTintUp;
uniform vec3 skyTintDown;

uniform float radianceClamp;
uniform float skipProbability;
uniform float shadowStrength;
uniform bool maxStepsIsMiss;
uniform int wavelengthSamples;

uniform float metalRoughness;
uniform vec3 metalSpecAlbedoRGB;
uniform float dieleRoughness;
uniform vec3 dieleAbsorptionRGB;
uniform vec3 dieleSpecAlbedoRGB;
uniform vec3 surfaceDiffuseAlbedoRGB;
uniform vec3 surfaceSpecAlbedoRGB;
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
#define RADIANCE_EPSILON 1.0e-6
#define M_PI 3.141592653589793

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2
#define MAT_VOLUM  3


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
    if (wiL.z < 0.0) return 0.0;
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

// ****************************        Dielectric        ****************************

#ifdef HAS_DIELECTRIC

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
    vec3 x = safe_normalize(cross(cross(n, wt), n));
    float sint = dot(wt, x);
    float sini = eta * sint; // Snell's law
    float sini2 = sini*sini;
    if (sini2 >= 1.0) return false; // Total internal reflection
    float cosi2 = (1.0 - sini*sini);
    float cosi = max(0.0, sqrt(cosi2));
    wi = -cosi*n + sini*x;
    return true;
}

RadianceType DIELECTRIC_SPEC_REFL_EVAL(in vec3 X, in vec3 woL, in Basis basis, in vec3 rgb)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 reflRGB = DIELECTRIC_SPECULAR_REFLECTANCE(dieleSpecAlbedoRGB, X, basis.nW, woW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType evaluateDielectric( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 rgb )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    RadianceType dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, rgb);
    RadianceType Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 h;
    float eta; // IOR ratio, et/ei
    if (reflected)
    { // Compute reflection half-vector
        h = safe_normalize(wiL + woL);
    }
    else
    { // Compute refraction half-vector
        bool entering = cosTheta(wiL) > 0.0;
        eta = entering ? ior : 1.0/ior;
        h = safe_normalize(wiL + eta*woL);
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    RadianceType f;
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
        f = (RadianceType(1.0) - Fr) * G * D * abs(im) * dwh_dwo / max(abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    }
    return f;
}

float pdfDielectric( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 rgb )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;
    RadianceType dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, rgb);
    RadianceType Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 h;
    float dwh_dwo;
    float pdf;
    if (reflected)
    {
        h = safe_normalize(wiL + woL);
        dwh_dwo = 1.0 / max(4.0*dot(woL, h), DENOM_TOLERANCE);
        pdf = averageComponent(Fr);
    }
    else
    { // Compute reflection half-vector
        bool entering = cosTheta(wiL) > 0.0;
        float eta = entering ? ior : 1.0/max(ior, DENOM_TOLERANCE);
        vec3 h = safe_normalize(wiL + eta*woL);
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        dwh_dwo = eta*eta * abs(om) / max(sqrtDenom*sqrtDenom, DENOM_TOLERANCE);
        pdf = 1.0 - averageComponent(Fr);
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    pdf *= microfacetPDF(h, roughness);
    return abs(pdf * dwh_dwo);
}

RadianceType sampleDielectric( in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 rgb,
                               inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_DIELE(wavelength_nm);
    RadianceType dielectricAlbedo = DIELECTRIC_SPEC_REFL_EVAL(X, woL, basis, rgb);
    RadianceType Fr = dielectricAlbedo * fresnelDielectricReflectance(woL.z, ior, 1.0);
    float roughness = DIELECTRIC_ROUGHNESS(dieleRoughness, X, basis.nW);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    float microPDF = microfacetPDF(m, roughness);
    float reflectProb = averageComponent(Fr);
    if (rand(rnd) < reflectProb)
    {
        // Compute specularly reflected ray direction
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute incident direction by reflecting woL about m
        if (wiL.z<0.0) wiL.z *= -1.0; // Reflect into positive hemisphere if necessary (ad hoc)
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        RadianceType f = Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
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
        if ( !refraction(ni, eta, woL, wiL) ) return RadianceType(0.0); // total internal reflection occurred
        wiL = -wiL; // As refract() computes the incident beam direction, and wiL is defined to be opposite to that.
        float cosi = dot(wiL, m);
        RadianceType Frm = RadianceType(fresnelDielectricReflectance(cosi, ior, 1.0));
        RadianceType Trm = RadianceType(1.0) - dielectricAlbedo * Frm;  // Fresnel transmittance
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
        RadianceType f = RadianceType(abs(im) * dwh_dwo * Trm * G * D) / max(abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
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

RadianceType METAL_SPEC_REFL_EVAL(in vec3 X, in vec3 woL, in Basis basis, in vec3 rgb)
{
    vec3 woW = localToWorld(woL, basis);
    vec3 reflRGB = METAL_SPECULAR_REFLECTANCE(metalSpecAlbedoRGB, X, basis.nW, woW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType evaluateMetal( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 rgb )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float Fr = fresnelMetalReflectance(woL.z, ior, k);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    RadianceType specAlbedo = METAL_SPEC_REFL_EVAL(X, woL, basis, rgb);
    RadianceType f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    return f;
}

float pdfMetal( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 rgb )
{
    float ior = IOR_DIELE(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    vec3 h = safe_normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float roughness = METAL_ROUGHNESS(metalRoughness, X, basis.nW);
    float pdf = microfacetPDF(h, roughness) * dwh_dwo;
    return pdf;
}

RadianceType sampleMetal( in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 rgb,
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
    RadianceType specAlbedo = METAL_SPEC_REFL_EVAL(X, woL, basis, rgb);
    RadianceType f = specAlbedo * Fr * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    float dwh_dwo; // Jacobian of the half-direction mapping
    dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
    pdfOut = microfacetPDF(m, roughness) * dwh_dwo;
    return f;
}

#endif

// ****************************        Surface        ****************************

#ifdef HAS_SURFACE

RadianceType SURFACE_DIFFUSE_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 woW, in vec3 rgb)
{
    vec3 reflRGB = SURFACE_DIFFUSE_REFLECTANCE(surfaceDiffuseAlbedoRGB, X, nW, woW);
    return rgbToAlbedo(reflRGB, rgb);
}

RadianceType SURFACE_SPEC_REFL_EVAL(in vec3 X, in vec3 nW, in vec3 woW, in vec3 rgb)
{
    vec3 reflRGB = SURFACE_SPECULAR_REFLECTANCE(surfaceSpecAlbedoRGB, X, nW, woW);
    return rgbToAlbedo(reflRGB, rgb);
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

RadianceType evaluateSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 rgb)
{
    vec3 woW = localToWorld(woL, basis);
    RadianceType diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, rgb);
    RadianceType    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, rgb);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    RadianceType f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    float E = fresnelDielectricReflectanceFast(woL.z, ior);
    f += (1.0 - E) * diffuseAlbedo/M_PI;
    return f;
}

float pdfSurface(in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in float wavelength_nm, in vec3 rgb)
{
    vec3 woW = localToWorld(woL, basis);
    RadianceType diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, rgb);
    RadianceType    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, rgb);
    float ior = surfaceIor;
    float E = fresnelDielectricReflectanceFast(woL.z, ior);
    float specWeight    = (E      )*averageComponent(specAlbedo);
    float diffuseWeight = (1.0 - E)*averageComponent(diffuseAlbedo);
    float weightSum = max(specWeight + diffuseWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float diffusePdf = pdfHemisphere(wiL);
    vec3 h = safe_normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float specularPdf = microfacetPDF(h, roughness) * dwh_dwo;
    return specProb*specularPdf + (1.0-specProb)*diffusePdf;
}

RadianceType sampleSurface(in vec3 X, in Basis basis, in vec3 woL, in float wavelength_nm, in vec3 rgb,
                           inout vec3 wiL, inout float pdfOut, inout vec4 rnd)
{
    vec3 woW = localToWorld(woL, basis);
    RadianceType diffuseAlbedo = SURFACE_DIFFUSE_REFL_EVAL(X, basis.nW, woW, rgb);
    RadianceType    specAlbedo = SURFACE_SPEC_REFL_EVAL(X, basis.nW, woW, rgb);
    float ior = surfaceIor;
    float roughness = SURFACE_ROUGHNESS(surfaceRoughness, X, basis.nW);
    float E = fresnelDielectricReflectanceFast(woL.z, ior);
    float specWeight    = (E      )*averageComponent(specAlbedo);
    float diffuseWeight = (1.0 - E)*averageComponent(diffuseAlbedo);
    float weightSum = max(specWeight + diffuseWeight, DENOM_TOLERANCE);
    float specProb = specWeight/weightSum;
    if (rand(rnd) >= specProb) // diffuse term, happens with probability 1-specProb
    {
        wiL = sampleHemisphere(rnd, pdfOut);
        float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
        pdfOut *= (1.0-specProb);
        return (1.0 - E) * diffuseAlbedo/M_PI;
    }
    else
    {
        vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
        pdfOut = specProb;
        if (wiL.z<DENOM_TOLERANCE) return RadianceType(0.0);
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float Fr = fresnelDielectricReflectanceFast(wiL.z, ior);
        RadianceType f = Fr * specAlbedo * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut *= microfacetPDF(m, roughness) * dwh_dwo;
        return f;
    }
}

#endif

// ****************************        BSDF common interface        ****************************

RadianceType evaluateBsdf( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 rgb,
                           inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    evaluateSurface(X, basis, woL, wiL, wavelength_nm, rgb); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      evaluateMetal(X, basis, woL, wiL, wavelength_nm, rgb); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return evaluateDielectric(X, basis, woL, wiL, wavelength_nm, rgb); }
#endif
}

RadianceType sampleBsdf( in vec3 X, in Basis basis, in vec3 woL, in int material, in float wavelength_nm, in vec3 rgb,
                         inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    sampleSurface(X, basis, woL, wavelength_nm, rgb, wiL, pdfOut, rnd); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      sampleMetal(X, basis, woL, wavelength_nm, rgb, wiL, pdfOut, rnd); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return sampleDielectric(X, basis, woL, wavelength_nm, rgb, wiL, pdfOut, rnd); }
#endif
}

float pdfBsdf( in vec3 X, in Basis basis, in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 rgb )
{
#ifdef HAS_SURFACE
    if (material==MAT_SURFA) { return    pdfSurface(X, basis, woL, wiL, wavelength_nm, rgb); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { return      pdfMetal(X, basis, woL, wiL, wavelength_nm, rgb); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { return pdfDielectric(X, basis, woL, wiL, wavelength_nm, rgb); }
#endif
}

#endif // HAS_GEOMETRY

#ifdef HAS_VOLUME

RadianceType VOLUME_ALBEDO_EVAL(in vec3 X, in vec3 rgb)
{
    vec3 sigma_s = VOLUME_SCATTERING_COLOR(volumeScatteringColorRGB, X);
    vec3 sigma_a = VOLUME_ABSORPTION_COLOR(volumeAbsorptionColorRGB, X);
    vec3 albedo_rgb = sigma_s / max(sigma_s + sigma_a, DENOM_TOLERANCE);
    return rgbToAlbedo(albedo_rgb, rgb);
}

RadianceType VOLUME_EXTINCTION_EVAL(in vec3 X, in vec3 rgb)
{
    vec3 sigma_s = VOLUME_SCATTERING_COLOR(volumeScatteringColorRGB, X);
    vec3 sigma_a = VOLUME_ABSORPTION_COLOR(volumeAbsorptionColorRGB, X);
    RadianceType extinction = volumeExtinction * RadianceType(VOLUME_EXTINCTION(volumeExtinction, X));
    vec3 sigma_t = extinction * (sigma_s + sigma_a) / lengthScale;
    return rgbToAlbedo(sigma_t, rgb);
}

RadianceType VOLUME_EXTINCTION_MAX_EVAL()
{
    return RadianceType(VOLUME_EXTINCTION_MAX(volumeExtinction));
}

float VOLUME_ANISOTROPY_EVAL(in vec3 X)
{
    return VOLUME_ANISOTROPY(volumeAnisotropy, X);
}

#endif // HAS_VOLUME

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
    if (material==MAT_SURFA) { N = vec3(   SDF_SURFACE(xyyp) -    SDF_SURFACE(xyyn),    SDF_SURFACE(yxyp) -    SDF_SURFACE(yxyn),    SDF_SURFACE(yyxp) -    SDF_SURFACE(yyxn)); return safe_normalize(N); }
#endif
#ifdef HAS_METAL
    if (material==MAT_METAL) { N = vec3(     SDF_METAL(xyyp) -      SDF_METAL(xyyn),      SDF_METAL(yxyp) -      SDF_METAL(yxyn),      SDF_METAL(yyxp) -      SDF_METAL(yyxn)); return safe_normalize(N); }
#endif
#ifdef HAS_DIELECTRIC
    if (material==MAT_DIELE) { N = vec3(SDF_DIELECTRIC(xyyp) - SDF_DIELECTRIC(xyyn), SDF_DIELECTRIC(yxyp) - SDF_DIELECTRIC(yxyn), SDF_DIELECTRIC(yyxp) - SDF_DIELECTRIC(yyxn)); return safe_normalize(N); }
#endif
#ifdef HAS_VOLUME
    if (material==MAT_VOLUM) { N = vec3(    SDF_VOLUME(xyyp) -     SDF_VOLUME(xyyn),     SDF_VOLUME(yxyp) -     SDF_VOLUME(yxyn),     SDF_VOLUME(yyxp) -     SDF_VOLUME(yyxn)); return safe_normalize(N); }
#endif
}

#ifdef HAS_NORMALMAP
void perturbNormal(in vec3 X, in Basis basis, int material, inout vec3 nW)
{
#ifdef HAS_SURFACE_NORMALMAP
    if (material==MAT_SURFA) 
    {
        vec3 nL = SURFACE_NORMAL_MAP(X);
        nW = localToWorld(normalize(2.0*nL - vec3(1.0)), basis);
    }
#endif
#ifdef HAS_METAL_NORMALMAP
    if (material==MAT_METAL) 
    {
        vec3 nL = METAL_NORMAL_MAP(X);
        nW = localToWorld(normalize(2.0*nL - vec3(1.0)), basis);
    }
#endif
#ifdef HAS_DIELECTRIC_NORMALMAP
    if (material==MAT_DIELE) 
    {
        vec3 nL = DIELECTRIC_NORMAL_MAP(X);
        nW = localToWorld(normalize(2.0*nL - vec3(1.0)), basis);
    }
#endif
}
#endif

// MC-estimates the amount of light transmitted along an infinite ray
RadianceType Visibility(in vec3 pW, in vec3 rayDir, in vec3 rgb, inout vec4 rnd, bool inVolume)
{
    vec3 pW_next;
    int hitMaterial;
    bool hit = traceRay(pW, rayDir, pW_next, hitMaterial, maxLengthScale);
#ifndef HAS_VOLUME
    return RadianceType(!hit);
#else

    // Deal with all volume interactions until the ray either hits a surface or escapes
    const int maxVolHits = 16;
    int volHits = 0;
    float eps = 3.0*minLengthScale;
    RadianceType sigma_t_max = VOLUME_EXTINCTION_MAX_EVAL();
    RadianceType inv_sigma_t_max = 1.0/max(sigma_t_max, DENOM_TOLERANCE);
#ifndef DISPERSION_ENABLED
    // Choose color channel to sample
    int channel = randomChannel(rnd);
    float inv_sigma_t_max_c = getChannel(inv_sigma_t_max, channel);
#endif
    while ((inVolume || hitMaterial==MAT_VOLUM) && volHits<maxVolHits)
    {
        // If ray lies in the volume interior
        if (inVolume)
        {
            // Do Woodcock tracking to check for a possible volume scattering event along the in-volume segment pW -> pW_next
            // (except if the current ray lies inside a dielectric, in which case there is no volumetric scattering, i.e. the dielectric displaces the volume)
            vec3 pScatter;
            RadianceType sigma_t;
            bool interacted = false;
            float segmentLength = length(pW_next - pW);
#ifdef DISPERSION_ENABLED
            float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max;
#else
            float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max_c;           
#endif

            int steps = 0;
            while (distanceMarched<segmentLength-eps && steps<__MAX_VOLUME_STEPS__)
            {
                pScatter = pW + distanceMarched*rayDir;
                sigma_t = VOLUME_EXTINCTION_EVAL(pScatter, rgb);
#ifdef DISPERSION_ENABLED
                interacted = bool(rand(rnd) < sigma_t * inv_sigma_t_max);
                if (interacted) break;
                distanceMarched += -log(rand(rnd)) * inv_sigma_t_max;
#else
                float sigma_t_c = getChannel(sigma_t, channel);
                interacted = bool(rand(rnd) < sigma_t_c * inv_sigma_t_max_c);
                if (interacted) break;
                distanceMarched += -log(rand(rnd)) * inv_sigma_t_max_c;
#endif
                steps++;
            }

            if (interacted) return RadianceType(0.0); // on absorption/scattering, return zero visibility
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
                    return RadianceType(0.0);
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
    bool visible = (!hit || volHits==maxVolHits);
#ifdef DISPERSION_ENABLED
    return float(visible); // if no hit or scattering, return visibility 1.0 (transparent)
#else
    vec3 V = vec3(0.0);
    setChannel(V, channel, float(visible));
    return V;
#endif

#endif // HAS_VOLUME
}


////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

vec3 environmentRadianceRGB(in vec3 dir)
{
    float rot_phi = M_PI*envMapPhiRotation/180.0;
    float rot_theta = M_PI*envMapThetaRotation/180.0;
    float phi = atan(dir.x, dir.z) + M_PI + rot_phi;
    phi -= 2.0*M_PI*floor(phi/(2.0*M_PI)); // wrap phi to [0, 2*pi]
    float theta = mod(acos(dir.y) + rot_theta, M_PI);             // theta in [0, pi]
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 RGB = vec3(1.0);
    if (haveEnvMap)
    {
        RGB = texture(envMap, vec2(u,v)).rgb;
    }
    float St = sin(rot_theta);
    vec3 color_pole = vec3(cos(rot_phi)*St, cos(rot_theta), sin(rot_phi)*St);
    float t = dot(dir, color_pole);
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

vec3 sampleSunDir(inout vec4 rnd, inout float pdf)
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
    pdf = 1.0/solid_angle;
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

float sunPowerEstimate(Basis basis)
{
    float theta_max = sunAngularSize * M_PI/180.0;
    float vis_factor = clamp(dot(basis.nW, sunDir) + sin(theta_max), 0.0, 1.0);
    float solid_angle = 2.0*M_PI*(1.0 - cos(theta_max));
    return vis_factor * solid_angle * sunPower * maxComponent(sunColor);
}

#ifdef HAS_GEOMETRY

RadianceType sampleLightAtSurface(Basis basis, in vec3 rgb, inout vec4 rnd, inout vec3 wiL, inout vec3 wiW, inout float lightPdf)
{
    if (sunPower<RADIANCE_EPSILON && skyPower<RADIANCE_EPSILON) 
        return RadianceType(0.0);
    
    // Light sampling (choose either sun or sky)
    float sunWeight = sunPower;
    float skyWeight = skyPower;
    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    float skyProb = max(PDF_EPSILON, 1.0-sunProb);
    bool chooseSun = (rand(rnd) <= sunProb);
    float skyPdf, sunPdf;
    RadianceType Li = RadianceType(0.0);
    if (chooseSun)
    {
        wiW = sampleSunDir(rnd, sunPdf);
        wiL = worldToLocal(wiW, basis);
        lightPdf = sunProb*sunPdf;
        if (wiL.z < 0.0) return RadianceType(0.0);
        Li += sunRadiance(wiW, rgb);
    }
    else
    {
        wiL = sampleHemisphere(rnd, skyPdf);
        wiW = localToWorld(wiL, basis);
        lightPdf = skyProb*skyPdf;
        Li += environmentRadiance(wiW, rgb);
    }
    return Li;
}

// Estimate direct radiance at the given surface vertex
RadianceType directSurfaceLighting(in vec3 pW, Basis basis, in vec3 woW, in int material,
                                   float wavelength_nm, in vec3 rgb, inout vec4 rnd, bool inVolume, inout float lightPdf)
{
    vec3 wiL, wiW; // direction of sampled direct light (*towards* the light)
    RadianceType Li = sampleLightAtSurface(basis, rgb, rnd, wiL, wiW, lightPdf);
    if ( averageComponent(Li) < RADIANCE_EPSILON ) return RadianceType(0.0);
    vec3 dPw = 3.0*minLengthScale * basis.nW;
    RadianceType V = Visibility(pW+dPw, wiW, rgb, rnd, inVolume);
    Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-V));
    // Apply MIS weight with the BSDF pdf for the sampled direction
    vec3 woL = worldToLocal(woW, basis);
    float bsdfPdf = pdfBsdf(pW, basis, woL, wiL, material, wavelength_nm, rgb);
    if (bsdfPdf<PDF_EPSILON) return RadianceType(0.0);
    RadianceType f = evaluateBsdf(pW, basis, woL, wiL, material, wavelength_nm, rgb, rnd);
    float misWeight = powerHeuristic(lightPdf, bsdfPdf);
    RadianceType fOverPdf = min(RadianceType(radianceClamp), f/max(PDF_EPSILON, lightPdf));
    return fOverPdf * Li * abs(dot(wiW, basis.nW)) * misWeight;
}

#endif // HAS_GEOMETRY

#ifdef HAS_VOLUME

RadianceType sampleLightInVolume(in vec3 rgb, inout vec4 rnd, inout vec3 wiW, inout float lightPdf)
{
    // Light sampling (choose either sun or sky)
    float sunWeight = sunPower;
    float skyWeight = skyPower;
    float sunProb = clamp(max(sunWeight, DENOM_TOLERANCE) / max(skyWeight+sunWeight, DENOM_TOLERANCE), 0.0, 1.0);
    float skyProb = max(PDF_EPSILON, 1.0-sunProb);
    bool chooseSun = (rand(rnd) <= sunProb);
    float skyPdf, sunPdf;
    RadianceType Li = RadianceType(0.0);
    if (chooseSun)
    {
        wiW = sampleSunDir(rnd, sunPdf);
        lightPdf = sunProb*sunPdf;
        Li += sunRadiance(wiW, rgb);
    }
    else
    {
        wiW = sampleSphere(rnd, skyPdf);
        lightPdf = skyProb*skyPdf;
        Li += environmentRadiance(wiW, rgb);
    }
    return Li;
}

// Estimate direct radiance at the given volumetric vertex
RadianceType directVolumeLighting(in vec3 pW, in vec3 woW, in vec3 rgb, inout vec4 rnd, bool inVolume)
{
    vec3 wiW; // direction of sampled direct light (*towards* the light)
    float lightPdf;
    RadianceType Li = sampleLightInVolume(rgb, rnd, wiW, lightPdf);
    RadianceType V = Visibility(pW, wiW, rgb, rnd, inVolume);
    Li *= abs(RadianceType(1.0) - shadowStrength*(RadianceType(1.0)-V));
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

RadianceType samplePath(in vec3 primaryStart, in vec3 primaryDir, 
                        float wavelength_nm, in vec3 rgb, inout vec4 rnd)

{
    // Perform pathtrace to estimate the primary ray radiance, L
    RadianceType L = RadianceType(0.0);
    float misWeight = 1.0; // For MIS book-keeping
    vec3 pW = primaryStart;
    vec3 rayDir = primaryDir; // (opposite to light direction)
    bool inVolume = false;
#ifdef HAS_VOLUME
    RadianceType sigma_t_max = VOLUME_EXTINCTION_MAX_EVAL();
    RadianceType inv_sigma_t_max = 1.0/sigma_t_max;
    inVolume = SDF_VOLUME(primaryStart) < 0.0;
#endif
    float eps = 3.0*minLengthScale;
#ifdef HAS_DIELECTRIC
    bool inDielectric = SDF_DIELECTRIC(primaryStart) < 0.0;
#endif

    RadianceType throughput = RadianceType(1.0);

#ifndef DISPERSION_ENABLED
    // choose volume color channel to sample for this entire path
    int volume_channel = randomChannel(rnd);
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
#ifndef DISPERSION_ENABLED
                // channel to sample
                float inv_sigma_t_max_c = getChannel(inv_sigma_t_max, volume_channel);
#endif
                // Do Woodcock tracking to check for a possible volume absorption/scattering event along the in-volume segment pW -> pW_next
                // (except if the current ray lies inside a dielectric, in which case there is no volumetric scattering, i.e. the dielectric displaces the volume)
                vec3 pScatter;
                RadianceType sigma_t;
                bool interacted = false;
#ifdef HAS_DIELECTRIC
                if (!inDielectric)
#endif
                {
#ifdef DISPERSION_ENABLED
                    float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max;
#else
                    float distanceMarched = -log(rand(rnd)) * inv_sigma_t_max_c;
#endif
                    int steps = 0;
                    while (distanceMarched<segmentLength-eps && steps<__MAX_VOLUME_STEPS__)
                    {
                        pScatter = pW + distanceMarched*rayDir;
                        sigma_t = VOLUME_EXTINCTION_EVAL(pScatter, rgb);
#ifdef DISPERSION_ENABLED
                        interacted = bool(rand(rnd) < sigma_t * inv_sigma_t_max);
                        if (interacted) break;
                        distanceMarched += -log(rand(rnd)) * inv_sigma_t_max;
#else
                        float sigma_t_c = getChannel(sigma_t, volume_channel);
                        interacted = bool(rand(rnd) < sigma_t_c * inv_sigma_t_max_c);
                        if (interacted) break;
                        distanceMarched += -log(rand(rnd)) * inv_sigma_t_max_c;
#endif
                        steps++;
                    }
                }

                // Absorption/scattering event
                if (interacted)
                {
                    RadianceType albedo = VOLUME_ALBEDO_EVAL(pScatter, rgb);
#ifdef DISPERSION_ENABLED
                    absorb = (rand(rnd) > albedo);
#else
                    float albedo_c = getChannel(albedo, volume_channel);
                    absorb = (rand(rnd) > albedo_c);
#endif
                    if (absorb)
                    {
#ifdef HAS_VOLUME_EMISSION
                        RadianceType emission = VOLUME_EMISSION_EVAL(pScatter, rgb);
                        L += throughput * emission; // add volume emission term
#endif
#ifndef DISPERSION_ENABLED
                        setChannel(throughput, volume_channel, 0.0);
#endif
                        break;
                    }

                    L += throughput * directVolumeLighting(pScatter, -rayDir, rgb, rnd, inVolume); // add contribution due to direct lighting at vertex
                    throughput *= albedo; // update throughput due to scattering
                    if (maxComponent(throughput) < THROUGHPUT_EPSILON) break;
                    pW = pScatter; // continue to next vertex
                    rayDir = samplePhaseFunction(rayDir, rnd); // sample scattered dir (NB, PF is the angle PDF so the MC PDF denom. cancels it in the L estimator)
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

        if (maxComponent(throughput) < THROUGHPUT_EPSILON) break;

        if (volHits == __MAX_SCATTERS__) hit = false;
#endif // HAS_VOLUME

        // If this ray didn't scatter and missed all geometry, add environment light term and terminate path
        if (!hit)
        {
            RadianceType Li = RadianceType(0.0);
            if (!(vertex==0 && !envMapVisible))      Li += environmentRadiance(rayDir, rgb);
            if (!(vertex==0 && !sunVisibleDirectly)) Li += sunRadiance(rayDir, rgb);
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
            RadianceType absorption = rgbToAlbedo(dieleAbsorptionRGB, rgb);
            throughput *= exp(-rayLength*absorption);
        }
#endif

        // This ray didn't scatter but hit some geometry, so deal with the surface interaction.
        // First, compute normal at current surface vertex
        pW = pW_next;
        vec3 nW = normal(pW, hitMaterial);
        Basis basis = makeBasis(nW);
#ifdef HAS_NORMALMAP
        perturbNormal(pW, basis, hitMaterial, nW);
        basis = makeBasis(nW);
#endif

        // Sample BSDF for the next bounce direction
        vec3 woW = -rayDir;
        vec3 woL = worldToLocal(woW, basis);
        vec3 wiL;
        float bsdfPdf;
        RadianceType f = sampleBsdf(pW, basis, woL, hitMaterial, wavelength_nm, rgb, wiL, bsdfPdf, rnd);
        vec3 wiW = localToWorld(wiL, basis);
        rayDir = wiW; // Update ray direction

        // Detect dielectric transmission
#ifdef HAS_DIELECTRIC
        if (hitMaterial==MAT_DIELE && dot(woW, nW)*dot(wiW, nW) < 0.0)
        {
            inDielectric = !inDielectric;
        }
#endif

#ifdef HAS_VOLUME_EMISSION
        // Add volumetric emission at the surface point, if present (treating it as isotropic radiance field)
        L += throughput * VOLUME_EMISSION_EVAL(pW, rgb);
#endif

        // Update path throughput
        RadianceType fOverPdf = min(RadianceType(radianceClamp), f/max(PDF_EPSILON, bsdfPdf));
        throughput *= fOverPdf * abs(dot(wiW, nW));
        if (maxComponent(throughput) < THROUGHPUT_EPSILON) break;

        // Add direct lighting term at current surface vertex
        float lightPdf = 0.0;
#ifdef HAS_DIELECTRIC
        if (!inDielectric)
#endif
        {
            L += throughput * directSurfaceLighting(pW, basis, woW, hitMaterial, wavelength_nm, rgb, rnd, inVolume, lightPdf);
        }

        // Prepare for tracing the bounce ray
        pW += nW * sign(dot(wiW, nW)) * eps;           // perturb vertex into half-space of scattered ray
        misWeight = powerHeuristic(bsdfPdf, lightPdf); // compute MIS weight for bounce ray
#endif // HAS_GEOMETRY
    }
    return L;
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
        // Jitter over pixel
        vec2 pixelj = pixel + (-0.5 + vec2(rand(rnd), rand(rnd)));

        // Compute world ray direction for this fragment
        vec3 primaryStart, primaryDir;
        constructPrimaryRay(pixelj, rnd, primaryStart, primaryDir);

        // Perform pathtrace to estimate the primary ray radiance, L
        L += samplePath(primaryStart, primaryDir, wavelength_nm, rgb, rnd);
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
