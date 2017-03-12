
uniform sampler2D PosData;
uniform sampler2D DirData;
uniform sampler2D RngData;
uniform sampler2D RgbData;

varying vec2 vTexCoord;
uniform float SceneScale;
uniform float roughness;


/////////////////////////////////////////////////////////////////////////
// Basis transforms
/////////////////////////////////////////////////////////////////////////

float CosTheta2(in vec3 nLocal) { return nLocal.z*nLocal.z; }
float TanTheta2(in vec3 nLocal)
{
	float ct2 = CosTheta2(nLocal);
	return max(0.0, 1.0 - ct2) / max(ct2, 1.0e-7);
}

float TanTheta(in vec3 nLocal)  { return sqrt(max(0.0, TanTheta2(nLocal))); }

void Basis(in vec3 nWorld, inout vec3 tWorld, inout vec3 bWorld)
{
	if (abs(nWorld.z) < abs(nWorld.x))
	{
		tWorld.x =  nWorld.z;
		tWorld.y =  0.0;
		tWorld.z = -nWorld.x;
	}
	else
	{
		tWorld.x =  0.0;
		tWorld.y = nWorld.z;
		tWorld.z = -nWorld.y;
	}
	tWorld = normalize(tWorld);
	bWorld = cross(nWorld, tWorld);
}

vec3 worldToLocal(in vec3 vWorld,
				  in vec3 nWorld, in vec3 tWorld, in vec3 bWorld)
{
	return vec3( dot(vWorld, tWorld),
			 	 dot(vWorld, bWorld),
				 dot(vWorld, nWorld) );
}

vec3 localToWorld(in vec3 vLocal,
				  in vec3 nWorld, in vec3 tWorld, in vec3 bWorld)
{
	return tWorld*vLocal.x + bWorld*vLocal.y + nWorld*vLocal.z;
}



/////////////////////////////////////////////////////////////////////////
// Beckmann Microfacet formulae
/////////////////////////////////////////////////////////////////////////

// m = the microfacet normal (in the local space where z = the macrosurface normal)
vec3 microfacetSample(inout vec4 rnd)
{
	float phiM = (2.0 * M_PI) * rand(rnd);
	float cosPhiM = cos(phiM);
	float sinPhiM = sin(phiM);
	float tanThetaMSqr = -roughness*roughness * log(max(1.0e-6, rand(rnd)));
	float cosThetaM = 1.0 / sqrt(1.0 + tanThetaMSqr);
	float sinThetaM = sqrt(max(0.0, 1.0 - cosThetaM*cosThetaM));
	return normalize(vec3(sinThetaM*cosPhiM, sinThetaM*sinPhiM, cosThetaM));
}

// Shadow-masking function
// Approximation from Walter et al (v = arbitrary direction, m = microfacet normal)
float smithG1(in vec3 vLocal, in vec3 mLocal)
{
	float tanTheta = abs(TanTheta(vLocal));
	if (tanTheta < 1.0e-6) 
	{
		return 1.0; // perpendicular incidence -- no shadowing/masking
	}
	// Back side is not visible from the front side, and the reverse.
	if (dot(vLocal, mLocal) * vLocal.z <= 0.0)
	{
		return 0.0;
	}
	// Rational approximation to the shadowing/masking function (Walter et al)  (<0.35% rel. error)
	float a = 1.0 / (roughness * tanTheta);
	if (a >= 1.6) 
	{
		return 1.0;
	}
	float aSqr = a*a;
	return (3.535*a + 2.181*aSqr) / (1.0 + 2.276*a + 2.577*aSqr);
}

float smithG2(in vec3 woL, in vec3 wiL, in vec3 mLocal)
{
	return smithG1(woL, mLocal) * smithG1(wiL, mLocal);
}


/////////////////////////////////////////////////////////////////////////
// Fresnel formulae
/////////////////////////////////////////////////////////////////////////


// D is direction of incident light
// N is normal pointing from medium to vacuum
float reflectionDielectric(in vec3 D, in vec3 N, float ior)
{
	float cosi = dot(-D, N);
	bool entering = cosi > 0.0;
	float ei, et;
	if (entering)
	{
		ei = 1.0; // Incident from vacuum, if entering
		et = ior; // Transmitted to internal medium, if entering
	}
	else
	{
		ei = ior;  // Incident from internal medium, if exiting
		et = 1.0;  // Transmitted to vacuum, if exiting
		N *= -1.0; // Flip normal (so normal is always opposite to incident light direction)
		cosi *= -1.0;
	}

	// Compute sint from Snells law:
	float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));

	// Handle total internal reflection
	if (sint >= 1.0) return 1.0;
	
	float cost = sqrt(max(0.0, 1.0 - sint*sint));
	float cosip = abs(cosi);
	const float epsilon = 1.0e-8;
	float rParallel      = (et*cosip - ei*cost) / max(et*cosip + ei*cost, epsilon);
	float rPerpendicular = (ei*cosip - et*cost) / max(ei*cosip + et*cost, epsilon);
	return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}


// D is direction of incident light
// N is normal pointing from medium to vacuum
// returns radiance multiplier (1 for reflection, varying for transmission)
float sampleDielectric(inout vec3 X, inout vec3 D, vec3 N, float ior, inout vec4 rnd)
{
	float normalEpsilon = 2.0e-5*SceneScale;
	float R = reflectionDielectric(D, N, ior);

	// Sample microfacet normal in z-local space
	vec3 mLocal = microfacetSample(rnd);

	vec3 T, B;
	Basis(N, T, B);
	vec3 mWorld = localToWorld(mLocal, N, T, B);

	float weight;

	// Reflection, with probability R
	if (R >= rand(rnd))
	{
		float cosi = dot(-D, mWorld);
		bool entering = cosi > 0.0;

		// displace ray start into reflected halfspace:
		if (entering)
		{
			X += normalEpsilon*N;
		}
		else
		{
			X -= normalEpsilon*N;
		}

		// Compute reflected direction by reflecting D about mWorld
		vec3 Dr = D - 2.0*mWorld*dot(D,mWorld);

		// Compute shadowing term
		vec3 DL  = worldToLocal(D, N, T, B);
		vec3 DrL = worldToLocal(Dr, N, T, B);
		float G = smithG2(DL, DrL, mLocal);

		// Compute Monte-Carlo sample weight
		// (From "Microfacet Models for Refraction through Rough Surfaces", Walter et. al, 2007)
		float dn = max(abs(dot(D,N)), 1.0e-6);
		float mn = max(abs(mLocal.z), 1.0e-6);
		weight = abs(cosi) * G  / (mn * dn);

		// Update direction
		D = Dr;
	}

	// Refraction, with probability (1-R)
	else 
	{
		float cosi = dot(-D, mWorld);
		bool entering = cosi > 0.0;

		// Compute transmitted (specularly refracted) ray direction
		float r;
		vec3 ni; // normal pointing into incident halfspace
		if (entering)
		{
			// Entering, incident halfspace is outside microfacet
			r = 1.0/max(ior, 1.0e-6);
			ni = mWorld;
			X -= normalEpsilon*N; // displace ray start into transmitted halfspace
		}
		else
		{
			// Exiting, incident halfspace is inside microfacet
			r = ior;
			ni = -mWorld;
			X += normalEpsilon*N; // displace ray start into transmitted halfspace
		}

		float eta = 1.0/max(r, 1.0e-6);

		// Compute sint from Snells law:
		float sint = r * sqrt(max(0.0, 1.0 - cosi*cosi));

		// sint<=1.0 guaranteed as total internal reflection already excluded
		float cost = sqrt(max(0.0, 1.0 - sint*sint)); 

		// Compute transmitted direction
		vec3 Dt = r*D + (r*abs(cosi) - cost)*ni; 

		// Compute shadowing term
		vec3 DL  = worldToLocal(D, N, T, B);
		vec3 DtL = worldToLocal(Dt, N, T, B);
		float G = smithG2(DL, DtL, mLocal);

		// Compute Monte-Carlo sample weight
		// (From "Microfacet Models for Refraction through Rough Surfaces", Walter et. al, 2007)
		float dn = max(abs(dot(D,N)), 1.0e-6);
		float mn = max(abs(mLocal.z), 1.0e-6);
		weight = abs(cosi) * G  / (mn * dn);

		// Update direction
		D = Dt;
	}

	return weight;
}


//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

SDF_FUNC

IOR_FUNC

SAMPLE_FUNC


//////////////////////////////////////////////////////////////
// SDF tracing
//////////////////////////////////////////////////////////////


vec3 NORMAL( in vec3 X )
{
	// Compute normal as gradient of SDF
	float normalEpsilon = 2.0e-5*SceneScale;
	vec3 eps = vec3(normalEpsilon, 0.0, 0.0);
	vec3 nor = vec3(
	    SDF_DIELE(X+eps.xyy) - SDF_DIELE(X-eps.xyy),
	    SDF_DIELE(X+eps.yxy) - SDF_DIELE(X-eps.yxy),
	    SDF_DIELE(X+eps.yyx) - SDF_DIELE(X-eps.yyx) );
	return normalize(nor);
}

void raytrace(inout vec4 rnd, 
			  inout vec3 X, inout vec3 D,
			  inout vec3 rgb, float wavelength)
{
	if (length(rgb) < 1.0e-8) return;
	float minMarchDist = 1.0e-5*SceneScale;
	float maxMarchDist = 1.0e5*SceneScale;
	float t = 0.0;
	float h = 1.0;
    for( int i=0; i<MAX_MARCH_STEPS; i++ )
    {
		if (h<minMarchDist || t>maxMarchDist) break;
		h = abs(SDF_DIELE(X + D*t));
		t += h;
    }
    X += t*D;
	if (t<maxMarchDist)
	{
		rgb *= SAMPLE(X, D, NORMAL(X), wavelength, rnd);
	}
	else
	{
		rgb *= 0.0; // terminate ray
	}
}


void main()
{
	vec3 X         = texture2D(PosData, vTexCoord).xyz;
	vec3 D         = texture2D(DirData, vTexCoord).xyz;
	vec4 rnd       = texture2D(RngData, vTexCoord);
	vec4 rgbw      = texture2D(RgbData, vTexCoord);

	float wavelength = 390.0 + (750.0 - 390.0)*rgbw.w;
	raytrace(rnd, X, D, rgbw.rgb, wavelength);

	float sgn = sign( SDF_DIELE(X) );

	gl_FragData[0] = vec4(X, sgn);
	gl_FragData[1] = vec4(D, 1.0);
	gl_FragData[2] = rnd;
	gl_FragData[3] = rgbw;
}


