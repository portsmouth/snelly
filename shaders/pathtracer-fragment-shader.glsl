
uniform sampler2D Radiance;         // 0
uniform sampler2D RngData;          // 1
uniform sampler2D WavelengthToRgb;  // 2
uniform sampler2D ICDF;             // 3

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
uniform float SceneScale;

uniform vec3 EmitterPos;
uniform vec3 EmitterDir;
uniform float EmitterRadius;
uniform float EmitterSpread; // in degrees

uniform float roughnessDiele;
uniform float roughnessMetal;

#define DENOM_TOLERANCE 1.0e-7

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

bool trace(in vec3 start, vec3 dir, inout vec3 hit, inout int material)
{
	float minMarchDist = 1.0e-5*SceneScale;
	float maxMarchDist = 1.0e5*SceneScale;
	float t = 0.0;
	float hDIELE = SceneScale;
	float hMETAL = SceneScale;
	float hDIFFU = SceneScale;
    for( int i=0; i<MAX_MARCH_STEPS; i++ )
    {
		if (hDIELE<minMarchDist) { material = MAT_DIELE; break; }
		if (hMETAL<minMarchDist) { material = MAT_METAL; break; }
		if (hDIFFU<minMarchDist) { material = MAT_DIFFU; break; }
		if (t>maxMarchDist) break;

		vec3 pW = start + t*dir;

		hDIELE = SDF_DIELE(pW);
		hMETAL = SDF_METAL(pW);
		hDIFFU = SDF_DIFFU(pW);
        t += min(hDIELE, min(hMETAL, hDIFFU));
    }
	if (t>maxMarchDist) return false;
	hit = start + t*dir;
	return true;
}

vec3 normal(in vec3 pW, int material)
{
	// Compute normal as gradient of SDF
	float normalEpsilon = 2.0e-5*SceneScale;
	vec3 e = vec3(normalEpsilon, 0.0, 0.0);
	vec3 xyyp = pW+e.xyy; vec3 xyyn = pW-e.xyy;
	vec3 yxyp = pW+e.yxy; vec3 yxyn = pW-e.yxy;
	vec3 yyxp = pW+e.yyx; vec3 yyxn = pW-e.yyx;
	vec3 N;
	if      (material==MAT_DIELE) { N = vec3(SDF_DIELE(xyyp)-SDF_DIELE(xyyn), SDF_DIELE(yxyp)-SDF_DIELE(yxyn), SDF_DIELE(yyxp) - SDF_DIELE(yyxn)); }
	else if (material==MAT_METAL) { N = vec3(SDF_METAL(xyyp)-SDF_METAL(xyyn), SDF_METAL(yxyp)-SDF_METAL(yxyn), SDF_METAL(yyxp) - SDF_METAL(yyxn)); }
	else                          { N = vec3(SDF_DIFFU(xyyp)-SDF_DIFFU(xyyn), SDF_DIFFU(yxyp)-SDF_DIFFU(yxyn), SDF_DIFFU(yyxp) - SDF_DIFFU(yyxn)); }
	return normalize(N);
}

bool Visible(in vec3 start, in vec3 end)
{
	float eps = 2.0e-5*SceneScale;
	vec3 dir = normalize(end - start);
	vec3 delta = eps * dir;
	start += delta;
	end   -= delta;
	vec3 p;
	int material;
	bool hit = trace(start, end, p, material);
	return !hit;
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

vec3 sampleHemisphere(inout vec4 rnd)
{
	// Do cosine-weighted sampling of hemisphere
	float r = sqrt(rand(rnd));
	float theta = 2.0 * M_PI * rand(rnd);
	float x = r * cos(theta);
	float y = r * sin(theta);
	float z = sqrt(1.0 - x*x - y*y);
	return vec3(x, y, z);
}


/////////////////////////////////////////////////////////////////////////
// Beckmann Microfacet formulae
/////////////////////////////////////////////////////////////////////////

// m = the microfacet normal (in the local space where z = the macrosurface normal)
float microfacetEval(in vec3 m, in float roughness)
{
	float tanTheta2 = tanTheta2(m);
	float cosTheta2 = cosTheta2(m);
	float roughnessSqr = roughness*roughness;
	float epsilon = 1.0e-5;
	float exponent = tanTheta2 / (roughnessSqr + epsilon);
	float D = exp(-exponent) / (M_PI * roughnessSqr * cosTheta2*cosTheta2);
	return D;
}

// m = the microfacet normal (in the local space where z = the macrosurface normal)
vec3 microfacetSample(inout vec4 rnd, in float roughness)
{
	float phiM = (2.0 * M_PI) * rand(rnd);
	float cosPhiM = cos(phiM);
	float sinPhiM = sin(phiM);
	float tanThetaMSqr = -roughness*roughness * log(max(1.0e-6, rand(rnd)));
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
	float a = 1.0 / (roughness * tanThetaAbs); // Rational approximation to the shadowing/masking function (Walter et al)  (<0.35% rel. error)
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
///		- dielectric <-> metal        (vacuum is a special case of dielectric)
///		- dielectric <-> diffuse      (vacuum is a special case of dielectric)
///		- dielectric <-> dielectric

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
		h = normalize(wiL + woL);
	}
	else
	{
		// Compute reflection half-vector
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

		// Microfacet distribution
		vec3 wh = m; // microfacet normal m = reflection half-vector
		float D = microfacetEval(wh, roughness);
		float G = smithG2(woL, wiL, m, roughness); // Shadow-masking function
		float f = R * D * G / (4.0*abs(cosTheta(wiL))*abs(cosTheta(woL)) + DENOM_TOLERANCE);

		// Return total BRDF and corresponding pdf
		pdfOut = microPDF * reflectProb;
		return f;
	}

	// transmission
	else
	{
		// Note, for transmission:
		//	woL.m < 0 means the transmitted light is entering the dielectric *from* the (vacuum) exterior of the microfacet
		//	woL.m > 0 means the transmitted light is exiting the dielectric *to* the (vacuum) exterior of the microfacet
		bool entering = dot(woL, m) < 0.0;

		// Compute transmitted (specularly refracted) ray direction
		float eta;
		vec3 ni; // normal pointing into incident halfspace
		if (entering)
		{
			// Entering, incident halfspace is outside microfacet
			eta = ior;
			ni = m;
		}
		else
		{
			// Exiting, incident halfspace is inside microfacet
			eta = 1.0/ior;
			ni = -m;
		}

		// Compute incident direction corresponding to known transmitted direction
		if ( !refraction(ni, eta, woL, wiL) )
		{
			return 0.0; // total internal reflection occurred
		}
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
	float tmp = (ior*ior + k*k) * cosi2;
	float twoEtaCosi = 2.0*ior*cosip;
	float Rparl2 = (tmp - twoEtaCosi + 1.0) / (tmp + twoEtaCosi + 1.0);
	float tmp_f = ior*ior + k*k;
	float Rperp2 = (tmp_f - twoEtaCosi + cosi2) / (tmp_f + twoEtaCosi + cosi2);
	return 0.5*(Rparl2 + Rperp2);
}

float evaluateMetal( in vec3 woL, in vec3 wiL, in float roughness, in float wavelength_nm )
{
	float ior = IOR_METAL(wavelength_nm);
	float k = K_METAL(wavelength_nm);
	return 0.0;
}

float pdfMetal( in vec3 woL, in vec3 wiL, in float roughness, in float wavelength_nm )
{
	float ior = IOR_METAL(wavelength_nm);
	float k = K_METAL(wavelength_nm);
	return 0.0;
}

float sampleMetal( in vec3 woL, in float roughness, in float wavelength_nm,
				   inout vec3 wiL, inout float pdfOut, inout vec4 rnd )
{
	float ior = IOR_METAL(wavelength_nm);
	float k = K_METAL(wavelength_nm);
	return 0.0;
}


// ****************************        Diffuse        ****************************

float evaluateDiffuse(in vec3 woL, in vec3 wiL)
{
	float reflectance = 0.5; // @todo: this will be in UI (like choice of dielectric & metal)
	return reflectance / M_PI;
}

float pdfDiffuse(in vec3 woL, in vec3 wiL)
{
	return abs(wiL.z)/M_PI;
}

float sampleDiffuse(in vec3 woL, inout vec3 wiL, 
				   inout float pdfOut, inout vec4 rnd)
{
	// Do cosine-weighted sampling of hemisphere
	wiL = sampleHemisphere(rnd);
	pdfOut = abs(wiL.z)/M_PI;
	float diffuseAlbedo = 0.5; // @todo: make this a wavelength-dependent color parameter
	return diffuseAlbedo / M_PI;
}

// ****************************        BSDF common interface        ****************************

float evaluateBsdf( in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm, 
	 			   inout vec4 rnd )
{
	return evaluateDielectric(woL, wiL, roughnessDiele, wavelength_nm);
	//if      (material==MAT_DIELE) { return evaluateDielectric(woL, wiL, roughnessDiele, wavelength_nm); }
	//else if (material==MAT_METAL) { return      evaluateMetal(woL, wiL, roughnessMetal, wavelength_nm); }
	//else                          { return    evaluateDiffuse(woL, wiL);                                }
}

float sampleBsdf( in vec3 woL, in int material, in float wavelength_nm,
				 inout vec3 wiL, inout float pdfOut, inout vec4 rnd ) 
{
	return sampleDielectric(woL, roughnessDiele, wavelength_nm, wiL, pdfOut, rnd);
	//if      (material==MAT_DIELE) { return sampleDielectric(woL, roughnessDiele, wavelength_nm, wiL, pdfOut, rnd); }
	//else if (material==MAT_METAL) { return      sampleMetal(woL, roughnessMetal, wavelength_nm, wiL, pdfOut, rnd); }
	//else                          { return    sampleDiffuse(woL,                                wiL, pdfOut, rnd); }
}

float pdfBsdf( in vec3 woL, in vec3 wiL, in int material, in float wavelength_nm )
{
	return pdfDielectric(woL, wiL, roughnessDiele, wavelength_nm);
	//if      (material==MAT_DIELE) { return pdfDielectric(woL, wiL, roughnessDiele, wavelength_nm); }
	//else if (material==MAT_METAL) { return      pdfMetal(woL, wiL, roughnessMetal, wavelength_nm); }
	//else                          { return    pdfDiffuse(woL, wiL);                                }
}


////////////////////////////////////////////////////////////////////////////////
// Pathtracing logic
////////////////////////////////////////////////////////////////////////////////


float environmentRadiance(in vec3 dir, float wavelength_nm)
{
	// For now, env-map is just a 'sky' which is a uniform blackbody with a specified temperature.
	const float sky_temp = 6500.0; // @todo: should be a param
	float boltzmann_factor = 1.43877737467e7 / (wavelength_nm*sky_temp);
    float l = wavelength_nm/360.0; // wavelength relative to 360nm
    return 1000.0 / (l*l*l*l*l*(exp(boltzmann_factor) - 1.0)) * cos(16.0*M_PI*dir.y) * 0.5*(1.0 + dir.x*dir.z);
}


bool emitterSample( in vec3 pW, in vec3 nW, 
			    	inout vec3 onLight, inout float lightPdf, inout vec4 rnd )
{	
	// project vertex onto emission plane
	float dPerp = dot(pW - EmitterPos, EmitterDir);
	vec3 pProj = pW - dPerp*EmitterDir;
	float dProj = length(pProj - EmitterPos);

	// thus radius of 'valid' disk containing possible source vertex on emission plane
	float spreadAngle = 0.5*abs(EmitterSpread)*M_PI/180.0;
	float rProjected = dPerp * tan(spreadAngle);
	if (dProj > rProjected + EmitterRadius) return false; // no direct light can reach the vertex

	// Choose a candidate point on the valid disk
	float rSample   = rProjected*sqrt(rand(rnd));
	float phiSample = 2.0*M_PI*rand(rnd);
	vec3 X = vec3(1.0, 0.0, 0.0);
	vec3 Z = vec3(0.0, 0.0, 1.0);
	vec3 u = cross(Z, EmitterDir);
	if ( length(u) < 1.0e-3 )
	{
		u = cross(X, EmitterDir);
	}
	u = normalize(u);
	vec3 v = cross(EmitterDir, u); // @todo, pass these u,v vectors to shader
	vec3 samplePos = rSample*(u*cos(phiSample) + v*sin(phiSample)); 
	if ( dot(samplePos, samplePos) > EmitterRadius*EmitterRadius )
	{
		// sample actually outside the emission disk
		return false;
	}

	onLight = EmitterPos + samplePos;

	// Compute solid-angle measure PDF of this sample (including samples which fell outside the disk)
	float diskArea = M_PI*rProjected*rProjected;
	vec3 toLight = onLight - pW;
	vec3 wiW = normalize(toLight);

	// No contribution if sampled light point is in opposite hemisphere to surface vertex:
	if (dot(wiW, nW) < 0.0) return false;

	float eps = 1.0e-6;
 	float jacobian = dot(toLight, toLight) / (abs(dot(EmitterDir, wiW)) + eps);
	lightPdf = jacobian / diskArea;
	return true;	
}

float powerHeuristic(const float a, const float b)
{
	float t = a*a;
	return t / (t + b*b);
}


float directLighting(in vec3 pW, Basis basis, in vec3 woW, in float emitterBrightness, 
	          		 in int material, float wavelength_nm, inout vec4 rnd)
{
	// choose whether to sample emitter or env map, (@todo: for now, just 50/50)
	float emissionProb = 0.0;

	float lightPdf;
	float Li;
	vec3 wiW;
	{
		// Env-map sampling
		//if ( rand(rnd) > emissionProb )
		{
			vec3 wiL = sampleHemisphere(rnd);
			lightPdf = (1.0 - emissionProb) * abs(wiL.z) / M_PI;
			wiW = localToWorld(wiL, basis);

			vec3 end = pW + wiW*SceneScale;
			bool visible = Visible(pW, end);
			if (!visible) Li = 0.0;
			else Li = environmentRadiance(wiW, wavelength_nm);
		}
		return Li;

		// emitter sampling
		/*
		else
		{
			// sample a point on the emission disk
			vec3 onLight;
			if ( !emitterSample(pW, basis.nW, onLight, lightPdf, rnd) ) Li = 0.0;
			else
			{
				// direction of direct light (*towards* the light)
				wiW = normalize(onLight - pW);

				// If light pointing away from vertex, or occluded, no direct light contribution.
				if (dot(wiW, basis.nW) < 0.0 || !Visible(pW, onLight)) Li = 0.0;
				else
				{
					lightPdf *= emissionProb;
					Li = emitterBrightness;
				}
			}
		}
		*/
	}
	
	// Apply MIS weight with the BSDF pdf for the sampled direction#
	vec3 woL = worldToLocal(woW, basis);
	vec3 wiL = worldToLocal(wiW, basis);

	float bsdfPdf = pdfBsdf(woL, wiL, material, wavelength_nm);
	const float PDF_EPSILON = 1.0e-5;
	if ( bsdfPdf<PDF_EPSILON ) return 0.0;

	float f = evaluateBsdf(woL, wiL, material, wavelength_nm, rnd);
	float misWeight = powerHeuristic(lightPdf, bsdfPdf);
	return f * Li * abs(dot(wiW, basis.nW)) * misWeight / max(PDF_EPSILON, lightPdf);
}


void main()
{
	vec4 rnd  = texture2D(RngData, vTexCoord);
	const float PDF_EPSILON = 1.0e-5;

	// Sample photon wavelength via the inverse CDF of the emission spectrum
	// (here w is spectral offset, i.e. wavelength = 360.0 + (750.0 - 360.0)*w)
	// (linear interpolation into the inverse CDF texture and RGB table should ensure smooth sampling over the range)
    float w = texture2D(ICDF, vec2(rand(rnd), 0.5)).r;
    float wavelength_nm = 360.0 + (750.0 - 360.0)*w;

    // @todo: we really need the physical 3-channel sensor sensitivies at the given wavelength
  	vec3 channelResponse = texture2D(WavelengthToRgb, vec2(w, 0.5)).rgb;

  	// Need something which gives the relative intensity of the emitter and the sky
  	float emitterBrightness = 1.0;

	// Jitter over pixel
	vec2 pixel = gl_FragCoord.xy;
	pixel += -0.5 + 0.5*vec2(rand(rnd), rand(rnd));

	// Compute world ray direction for this fragment
	vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
	float fh = camNear*tan(0.5*radians(camFovy)) / camZoom; // frustum height
	float fw = camAspect*fh;
	vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
	vec3 primaryDir = normalize(camNear*camDir + s); // ray direction

	// Raycast to first hit point
	vec3 pW;
	vec3 woW = -primaryDir;
	int material;
	bool hit = trace(camPos, primaryDir, pW, material);
	float zHit;

	float L;
	float throughput; 

	if ( !hit )
	{
		zHit = camFar;
		L = environmentRadiance(primaryDir, wavelength_nm);
		throughput = 1.0;
	}
	else
	{
		const int maxBounces = 4; // debug
		
		L = 0.0;
		throughput = 1.0;
		zHit = dot(pW - camPos, camDir);

		for (int bounce=0; bounce<maxBounces; ++bounce)
		{
			// Compute normal at current surface vertex
			vec3 nW = normal(pW, material);
			Basis basis = makeBasis(nW);

			// Add direct lighting term
			//L += throughput * directLighting(pW, basis, woW, emitterBrightness, material, wavelength_nm, rnd);

			// Sample BSDF for the next bounce direction
			vec3 woL = worldToLocal(woW, basis);
			vec3 wiL;
			float bsdfPdf;
			float f = sampleBsdf(woL, material, wavelength_nm, wiL, bsdfPdf, rnd);
			vec3 wiW = localToWorld(wiL, basis);

			// Update path throughput
			throughput *= f * abs(dot(wiW, nW)) / max(PDF_EPSILON, bsdfPdf);

			// Trace bounce ray
			float eps = 5.0e-5*SceneScale;
			pW += nW * dot(wiW, nW) * eps; // perturb vertex in direction of scattered ray
			vec3 pW_next;
			bool hit = trace(pW, wiW, pW_next, material);
			
			// Exit now if ray missed
			if (!hit)
			{
				float lightPdf = 1.0; //0.5;
				float misWeight = powerHeuristic(bsdfPdf, lightPdf);
				float Li = environmentRadiance(wiW, wavelength_nm);
				L += throughput * Li * misWeight;
				break;
			}

			// @todo: add BSDF-sampled emission term!
			{
				// check if ray segment (pW, pW_next) intersects the emitter disk.
				// if it does, add contribution (depending on spread) and terminate path.
			}

			// Update vertex
			vec3 rayDir = normalize(pW_next - pW);
			woW = -rayDir;
			pW = pW_next;
		}
	}

	vec3 color = channelResponse * L;

	float clipDepth = computeClipDepth(zHit, camNear, camFar);

	// Write updated radiance and sample count
	vec4 oldL = texture2D(Radiance, vTexCoord);
	float oldN = oldL.w;
	float newN = oldN + 1.0;
	vec3 newL = (oldN*oldL.rgb + color) / newN;

	gl_FragData[0] = vec4(newL, newN);
	gl_FragData[1] = rnd;
	gl_FragData[2] = pack_depth(clipDepth);
}








