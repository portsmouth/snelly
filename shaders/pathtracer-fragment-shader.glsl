
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
uniform int downRes;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camNear;
uniform float camFar;
uniform float camFovy; // degrees 
uniform float camZoom;
uniform float camAspect;
uniform float SceneScale;
uniform float SkyPower;
uniform vec3 diffuseAlbedoXYZ;
uniform float roughnessDiele;
uniform float roughnessMetal;
uniform vec3 absorptionDieleRGB;

#define DENOM_TOLERANCE 1.0e-7
#define HIT_TOLERANCE 1.0e-4
#define NORMAL_TOLERANCE 5.0e-4

#define MAT_VACUU  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_DIFFU  2

//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

SDF_FUNC

IOR_FUNC


///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit over specified segment
bool traceDistance(in vec3 start, in vec3 dir, float maxDist,
                   inout vec3 hit, inout int material)
{
    float minMarchDist = HIT_TOLERANCE*SceneScale;
    float sdf_diele = abs(SDF_DIELECTRIC(start));
    float sdf_metal = abs(SDF_METAL(start));
    float sdf_diffu = abs(SDF_DIFFUSE(start));
    float sdf = min(sdf_diele, min(sdf_metal, sdf_diffu));
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
        sdf_diffu = abs(SDF_DIFFUSE(pW));    if (sdf_diffu<minMarchDist) { material = MAT_DIFFU; break; }
        sdf = min(sdf_diele, min(sdf_metal, sdf_diffu));
        iters++;
    }
    hit = start + t*dir;
    if (t>=maxDist || iters>=MAX_MARCH_STEPS) return false;
    return true;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir, 
              inout vec3 hit, inout int material)
{
    float maxMarchDist = 2.0e2*SceneScale;
    material = MAT_VACUU;
    return traceDistance(start, dir, maxMarchDist, hit, material);
}


// (whether not occluded along finite length segment)
bool Visible(in vec3 start, in vec3 end)
{
    float eps = 20.0*HIT_TOLERANCE*SceneScale;
    vec3 dir = normalize(end - start);
    float maxDist = length(end - start);
    vec3 delta = eps * dir;
    vec3 hit;
    int material;
    bool occluded = traceDistance(start+delta, dir, maxDist, hit, material);
    return !occluded;
}

// (whether occluded along infinite ray)
bool Occluded(in vec3 start, in vec3 dir)
{
    float eps = 20.0*HIT_TOLERANCE*SceneScale;
    vec3 delta = eps * dir;
    vec3 p;
    int material;
    vec3 hit;
    bool occluded = traceRay(start+delta, dir, hit, material);
    return occluded;
}

vec3 normal(in vec3 pW, int material)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = NORMAL_TOLERANCE*SceneScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 xyyp = pW+e.xyy; vec3 xyyn = pW-e.xyy;
    vec3 yxyp = pW+e.yxy; vec3 yxyn = pW-e.yxy;
    vec3 yyxp = pW+e.yyx; vec3 yyxn = pW-e.yyx;
    vec3 N;
    if      (material==MAT_DIELE) { N = vec3(SDF_DIELECTRIC(xyyp)-SDF_DIELECTRIC(xyyn), SDF_DIELECTRIC(yxyp)-SDF_DIELECTRIC(yxyn), SDF_DIELECTRIC(yyxp) - SDF_DIELECTRIC(yyxn)); }
    else if (material==MAT_METAL) { N = vec3(SDF_METAL(xyyp)     -SDF_METAL(xyyn),      SDF_METAL(yxyp)     -SDF_METAL(yxyn),      SDF_METAL(yyxp)      - SDF_METAL(yyxn)); }
    else                          { N = vec3(SDF_DIFFUSE(xyyp)   -SDF_DIFFUSE(xyyn),    SDF_DIFFUSE(yxyp)   -SDF_DIFFUSE(yxyn),    SDF_DIFFUSE(yyxp)    - SDF_DIFFUSE(yyxn)); }
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
    vec3 nW;
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


/////////////////////////////////////////////////////////////////////////
// Sampling formulae
/////////////////////////////////////////////////////////////////////////

vec3 sampleHemisphere(inout vec4 rnd, inout float pdf)
{
    // Do cosine-weighted sampling of hemisphere
    float r = sqrt(rand(rnd));
    float theta = 2.0 * M_PI * rand(rnd);
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(1.0 - x*x - y*y);
    pdf = abs(z) / M_PI;
    return vec3(x, y, z);
}


float powerHeuristic(const float a, const float b)
{
    float t = a*a;
    return t / (t + b*b);
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
    return microfacetEval(m, roughness) * cosTheta(m);
}

// Shadow-masking function
// Approximation from Walter et al (v = arbitrary direction, m = microfacet normal)
float smithG1(in vec3 vLocal, in vec3 mLocal, float roughness)
{
    float tanThetaAbs = abs(tanTheta(vLocal));
    if (tanThetaAbs < 1.0e-6) return 1.0; // perpendicular incidence -- no shadowing/masking
    if (dot(vLocal, mLocal) * vLocal.z <= 0.0) return 0.0; // Back side is not visible from the front side, and the reverse.
    float epsilon = 1.0e-6;
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


/// @todo:
/// In general, have to deal with interfaces:
///     - dielectric <-> metal        (vacuum is a special case of dielectric)
///     - dielectric <-> diffuse      (vacuum is a special case of dielectric)
///     - dielectric <-> dielectric

/// Ray is either travelling in vacuum, or within the dielectric.


// ****************************        Dielectric        ****************************

/// Compute Fresnel reflectance at a dielectric interface (which has an "interior" and an "exterior").
/// cosi is the cosine to the (interior-to-exterior) normal of the incident ray
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

float evaluateDielectric( in vec3 woL, in vec3 wiL, in float roughness, in float wavelength_nm )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;

    // We call this when vertex is *on* the dielectric surface, which implies
    // in our case it is adjacent to vacuum. 
    float R = fresnelDielectricReflectance(woL.z, ior, 1.0);

    // Compute the reflection half-vector
    vec3 h;
    float eta; // IOR ratio, et/ei

    if (reflected)
    {
        // Compute reflection half-vector
        h = normalize(wiL + woL);
    }
    else
    {
        // Compute refraction half-vector
        bool entering = cosTheta(wiL) > 0.0;
        eta = entering ? ior : 1.0/ior;
        h = normalize(wiL + eta*woL);
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out

    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f;
    if (reflected)
    {
        f = R * D * G / (4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)) + DENOM_TOLERANCE);
    }
    else
    {
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        float dwh_dwo = eta*eta * abs(om) / (sqrtDenom*sqrtDenom);
        f = (1.0 - R) * G * D * abs(im) * dwh_dwo / (abs(cosTheta(wiL))*abs(cosTheta(woL)) + DENOM_TOLERANCE);
    }

    return f;
}

float pdfDielectric( in vec3 woL, in vec3 wiL, in float roughness, in float wavelength_nm )
{
    float ior = IOR_DIELE(wavelength_nm);
    bool reflected = cosTheta(wiL) * cosTheta(woL) > 0.0;

    // We call this when vertex is *on* the dielectric surface, which implies
    // in our case it is adjacent to vacuum. 
    float R = fresnelDielectricReflectance(woL.z, ior, 1.0);

    // Compute the reflection half-vector, and the Jacobian of the half-direction mapping
    vec3 h;
    float dwh_dwo;
    float pdf;

    if (reflected)
    {
        h = normalize(wiL + woL);
        dwh_dwo = 1.0 / (4.0*dot(woL, h) + DENOM_TOLERANCE);
        pdf = R;
    }
    else
    {
        // Compute reflection half-vector
        bool entering = cosTheta(wiL) > 0.0;
        float eta = entering ? ior : 1.0/ior;
        vec3 h = normalize(wiL + eta*woL);
        float im = dot(wiL, h);
        float om = dot(woL, h);
        float sqrtDenom = im + eta*om;
        dwh_dwo = eta*eta * abs(om) / (sqrtDenom*sqrtDenom + DENOM_TOLERANCE);
        pdf = 1.0 - R;
    }
    if (cosTheta(h)<0.0) h *= -1.0; // make sure half-vector points out

    pdf *= microfacetPDF(h, roughness);
    return abs(pdf * dwh_dwo);
}

float sampleDielectric( in vec3 woL, in float roughness, in float wavelength_nm,
                        inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_DIELE(wavelength_nm);

    // We call this when vertex is *on* the dielectric surface, which implies
    // in our case it is adjacent to vacuum. 
    float R = fresnelDielectricReflectance(woL.z, ior, 1.0);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    float microPDF = microfacetPDF(m, roughness);
    float reflectProb = R;

    // Choose whether to reflect or transmit randomly
    if (rand(rnd) < reflectProb)
    {
        // Compute specularly reflected ray direction
        wiL = -woL + 2.0*dot(woL, m)*m; // Compute incident direction by reflecting woL about m
        float D = microfacetEval(m, roughness);
        float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
        float f = R * D * G / (4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)) + DENOM_TOLERANCE);
        float dwh_dwo; // Jacobian of the half-direction mapping
        dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
        pdfOut = microPDF * reflectProb * dwh_dwo; // Return total BRDF and corresponding pdf
        return f;
    }

    // transmission
    else
    {
        // Note, for transmission:
        //  woL.m < 0 means the transmitted light is entering the dielectric *from* the (vacuum) exterior of the microfacet
        //  woL.m > 0 means the transmitted light is exiting the dielectric *to* the (vacuum) exterior of the microfacet
        bool entering = dot(woL, m) < 0.0;

        // Compute transmitted (specularly refracted) ray direction
        float eta;
        vec3 ni; // normal pointing into incident halfspace
        if (entering)
        {
            // Entering, incident halfspace is outside dielectric
            eta = ior; // = et/ei
            ni = m;
        }
        else
        {
            // Exiting, incident halfspace is inside dielectric
            eta = 1.0/ior; // = et/ei
            ni = -m;
        }

        // Compute incident direction corresponding to known transmitted direction
        if ( !refraction(ni, eta, woL, wiL) ) return 0.0; // total internal reflection occurred
        wiL = -wiL; // As refract() computes the incident beam direction, and wiL is defined to be opposite to that.
    
        // Compute Fresnel transmittance
        float cosi = dot(wiL, m);
        float R = fresnelDielectricReflectance(cosi, ior, 1.0);
        float T = 1.0 - R;

        // Evaluate microfacet distribution for the sampled half direction
        vec3 wh = m; // refraction half-vector = m
        float D = microfacetEval(wh, roughness);
        float G = smithG2(woL, wiL, wh, roughness); // Shadow-masking function
        float dwh_dwo; // Jacobian of the half-direction mapping
        float im = dot(wiL, m);
        {
            float om = dot(woL, m);
            float sqrtDenom = im + eta*om;
            dwh_dwo = eta*eta * abs(om) / (sqrtDenom*sqrtDenom + DENOM_TOLERANCE);
        }

        float f = abs(im) * dwh_dwo * T * G * D / (abs(cosTheta(wiL))*abs(cosTheta(woL)) + DENOM_TOLERANCE);
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

float evaluateMetal( in vec3 woL, in vec3 wiL, in float roughness, in float wavelength_nm )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float R = fresnelMetalReflectance(woL.z, ior, k);
    vec3 h = normalize(wiL + woL); // Compute the reflection half-vector
    float D = microfacetEval(h, roughness);
    float G = smithG2(woL, wiL, h, roughness);
    float f = R * D * G / max(4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)), DENOM_TOLERANCE);
    return f;
}

float pdfMetal( in vec3 woL, in vec3 wiL, in float roughness, in float wavelength_nm )
{
    float ior = IOR_DIELE(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float R = fresnelMetalReflectance(woL.z, ior, k);
    vec3 h = normalize(wiL + woL); // reflection half-vector
    float dwh_dwo = 1.0 / max(abs(4.0*dot(woL, h)), DENOM_TOLERANCE); // Jacobian of the half-direction mapping
    float pdf = microfacetPDF(h, roughness) * dwh_dwo;
    return pdf;
}

float sampleMetal( in vec3 woL, in float roughness, in float wavelength_nm,
                   inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
    float ior = IOR_METAL(wavelength_nm);
    float k = K_METAL(wavelength_nm);
    float R = fresnelMetalReflectance(woL.z, ior, k);
    vec3 m = microfacetSample(rnd, roughness); // Sample microfacet normal m
    wiL = -woL + 2.0*dot(woL, m)*m; // Compute wiL by reflecting woL about m
    float D = microfacetEval(m, roughness);
    float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
    float f = R * D * G / (4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)) + DENOM_TOLERANCE);
    float dwh_dwo; // Jacobian of the half-direction mapping
    dwh_dwo = 1.0 / max(abs(4.0*dot(woL, m)), DENOM_TOLERANCE);
    // Return total BRDF and corresponding pdf
    pdfOut = microfacetPDF(m, roughness) * dwh_dwo;
    return f;
}


// ****************************        Diffuse        ****************************

// @todo: albedo wavelength dependence?

float evaluateDiffuse(in vec3 woL, in vec3 wiL, in vec3 XYZ)
{
    vec3 albedo = diffuseAlbedoXYZ;
    return dot(albedo, XYZ) / M_PI;
}

float pdfDiffuse(in vec3 woL, in vec3 wiL)
{
    return abs(wiL.z) / M_PI;
}

float sampleDiffuse(in vec3 woL, in vec3 XYZ,
                    inout vec3 wiL, inout float pdfOut, inout vec4 rnd)
{
    // Do cosine-weighted sampling of hemisphere
    wiL = sampleHemisphere(rnd, pdfOut);
    vec3 albedo = diffuseAlbedoXYZ;
    return dot(albedo, XYZ) / M_PI;
}

// ****************************        BSDF common interface        ****************************

float evaluateBsdf( in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, in vec3 XYZ, 
                    inout vec4 rnd )
{
    if      (material==MAT_DIELE) { return evaluateDielectric(woL, wiL, roughnessDiele, wavelength_nm); }
    else if (material==MAT_METAL) { return      evaluateMetal(woL, wiL, roughnessMetal, wavelength_nm); }
    else                          { return    evaluateDiffuse(woL, wiL, XYZ);                                }
}

float sampleBsdf( in vec3 woL, in int material, in float wavelength_nm, in vec3 XYZ,
                  inout vec3 wiL, inout float pdfOut, inout vec4 rnd ) 
{
    if      (material==MAT_DIELE) { return sampleDielectric(woL, roughnessDiele, wavelength_nm, wiL, pdfOut, rnd); }
    else if (material==MAT_METAL) { return      sampleMetal(woL, roughnessMetal, wavelength_nm, wiL, pdfOut, rnd); }
    else                          { return    sampleDiffuse(woL,                 XYZ,           wiL, pdfOut, rnd); }
}

float pdfBsdf( in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm )
{
    if      (material==MAT_DIELE) { return pdfDielectric(woL, wiL, roughnessDiele, wavelength_nm); }
    else if (material==MAT_METAL) { return pdfMetal(woL, wiL, roughnessMetal, wavelength_nm);      }
    else                          { return pdfDiffuse(woL, wiL);                                   }
}


////////////////////////////////////////////////////////////////////////////////
// Light sampling
////////////////////////////////////////////////////////////////////////////////

float environmentRadiance(in vec3 dir, in vec3 XYZ)
{
    float phi = atan(dir.x, dir.z) + M_PI; // [0, 2*pi]
    float theta = acos(dir.y);             // [0, pi]
    float u = phi/(2.0*M_PI);
    float v = theta/M_PI;
    vec3 rgb = SkyPower * texture2D(envMap, vec2(u,v)).rgb; 

    // Here assuming texture is sRGB
    float X = 0.4124564*rgb.r + 0.3575761*rgb.g + 0.1804375*rgb.b;
    float Y = 0.2126729*rgb.r + 0.7151522*rgb.g + 0.0721750*rgb.b;
    float Z = 0.0193339*rgb.r + 0.1191920*rgb.g + 0.9503041*rgb.b;

    // convert to radiance via expansion in color matching functions
    vec3 c;
    c.x =  0.03382146*X - 0.02585410*Y - 0.00406490*Z;
    c.y = -0.02585410*X + 0.03209432*Y + 0.00227671*Z;
    c.z = -0.00406490*X + 0.00227671*Y + 0.00703345*Z;

    return dot(c, XYZ);
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
        bool occluded = Occluded(pW, wiW); 
        if (occluded) return 0.0;
        else 
            Li = environmentRadiance(wiW, XYZ);
    }

    // Apply MIS weight with the BSDF pdf for the sampled direction
    vec3 woL = worldToLocal(woW, basis);
    vec3 wiL = worldToLocal(wiW, basis);

    float bsdfPdf = pdfBsdf(woL, wiL, material, wavelength_nm);
    const float PDF_EPSILON = 1.0e-5;
    if ( bsdfPdf<PDF_EPSILON ) return 0.0;

    float f = evaluateBsdf(woL, wiL, material, wavelength_nm, XYZ, rnd);
    float misWeight = powerHeuristic(lightPdf, bsdfPdf);
    return f * Li * abs(dot(wiW, basis.nW)) * misWeight / max(PDF_EPSILON, lightPdf);
}



////////////////////////////////////////////////////////////////////////////////
// Pathtracing logic
////////////////////////////////////////////////////////////////////////////////

void pathtrace(vec2 pixel, vec4 rnd) // the current pixel
{
    const float PDF_EPSILON = 1.0e-5;

    // Sample photon wavelength via the inverse CDF of the emission spectrum
    // (here w is spectral offset, i.e. wavelength = 360.0 + (750.0 - 360.0)*w)
    // (linear interpolation into the inverse CDF texture and XYZ table should ensure smooth sampling over the range)
    float w = texture2D(ICDF, vec2(rand(rnd), 0.5)).r;
    float wavelength_nm = 390.0 + (750.0 - 390.0)*w;

    // Convert wavelength to XYZ tristimulus
    vec3 XYZ = texture2D(WavelengthToXYZ, vec2(w, 0.5)).rgb;

    // Jitter over pixel
    pixel += -0.5 + vec2(rand(rnd), rand(rnd));

    // Compute world ray direction for this fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
    float fh = camNear*tan(0.5*radians(camFovy)) / camZoom; // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    vec3 primaryDir = normalize(camNear*camDir + s); // ray direction

    // Raycast to first hit point
    vec3 pW;
    vec3 woW = -primaryDir;
    int rayMaterial = MAT_VACUU;
    int hitMaterial;
    bool hit = traceRay(camPos, primaryDir, pW, hitMaterial);
    float zHit;
    float L;
    float throughput; 


    if ( !hit )
    {
        zHit = camFar;
        L = environmentRadiance(primaryDir, XYZ);
        throughput = 1.0;
    }
    else
    {       
        L = 0.0;
        throughput = 1.0;
        zHit = dot(pW - camPos, camDir);

        for (int bounce=0; bounce<MAX_BOUNCES; ++bounce)
        {
            // Compute normal at current surface vertex
            vec3 nW = normal(pW, hitMaterial);
            Basis basis = makeBasis(nW);

            // Add direct lighting term
            L += throughput * directLighting(pW, basis, woW, hitMaterial, wavelength_nm, XYZ, rnd);

            // Sample BSDF for the next bounce direction
            vec3 woL = worldToLocal(woW, basis);
            vec3 wiL;
            float bsdfPdf;
            float f = sampleBsdf(woL, hitMaterial, wavelength_nm, XYZ, wiL, bsdfPdf, rnd);
            vec3 wiW = localToWorld(wiL, basis);

            // Update path throughput
            throughput *= f * abs(dot(wiW, nW)) / max(PDF_EPSILON, bsdfPdf);

            // Trace bounce ray
            float displacement = 20.0*HIT_TOLERANCE*SceneScale;
            float wiWnW = dot(wiW, nW);
            pW += nW * sign(wiWnW) * displacement; // perturb vertex into half-space of scattered ray
            vec3 pW_next;
            bool hit = traceRay(pW, wiW, pW_next, hitMaterial);

            // Exit now if ray missed
            if (!hit)
            {
                float hemispherePdf = abs(wiWnW);
                float lightPdf = hemispherePdf;
                float misWeight = powerHeuristic(bsdfPdf, lightPdf);
                float Li = environmentRadiance(wiW, XYZ);
                L += throughput * Li * misWeight;
                break;
            }

            // material in which ray propagates changes (only) on transmission
            vec3 incid_halfspace = dot(woW, nW) * nW;
            vec3 scatt_halfspace = dot(wiW, nW) * nW;
            if (dot(incid_halfspace, scatt_halfspace)<0.0) 
            {
                rayMaterial = hitMaterial;
            }

            // If the bounce ray lies inside a dielectric, apply Beer's law for absorption
            if (rayMaterial==MAT_DIELE)
            {
                float absorptionLength = length(pW_next - pW);
                throughput *= exp(-absorptionLength * dot(absorptionDieleRGB, RGB));
            }

            // Update vertex
            vec3 rayDir = normalize(pW_next - pW);
            woW = -rayDir;
            pW = pW_next;
        }
    }

    vec3 color = XYZ * L;
    //float clipDepth = computeClipDepth(zHit, camNear, camFar);

    // Write updated radiance and sample count
    vec4 oldL = texture2D(Radiance, vTexCoord);
    float oldN = oldL.w;
    float newN = oldN + 1.0;
    vec3 newL = (oldN*oldL.rgb + color) / newN;

    gl_FragData[0] = vec4(newL, newN);
    gl_FragData[1] = rnd;
    //gl_FragData[2] = pack_depth(clipDepth);
    
}

void RENDER_ALL()
{
    vec4 rnd = texture2D(RngData, vTexCoord);
    pathtrace(gl_FragCoord.xy, rnd);
}

void RENDER_BLOCKS()
{
    vec2 pixel = gl_FragCoord.xy - vec2(0.5, 0.5); // shift to get integer coords 
    float xrem = mod(float(pixel.x), float(downRes));
    float yrem = mod(float(pixel.y), float(downRes));
    float maxx = resolution.x-1.0;
    float maxy = resolution.y-1.0;
    vec4 rnd = texture2D(RngData, vTexCoord);

    // Sample only 1 pixel per block (the lower left one)
    float MOD_TOL = 1.0e-3;
    if (xrem<MOD_TOL && yrem<MOD_TOL)
    {
        // jitter pixel in block..
        pixel += float(downRes) * vec2(rand(rnd), rand(rnd));
        pixel.x = min(pixel.x, maxx);
        pixel.y = min(pixel.y, maxy);
        pathtrace(pixel, rnd);
    }
    else
    {
        gl_FragData[0] = texture2D(Radiance, vTexCoord);
        gl_FragData[1] = rnd;
        //gl_FragData[2] = pack_depth(clipDepth);
    }
}











