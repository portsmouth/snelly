
uniform sampler2D PosData;
uniform sampler2D DirData;
uniform sampler2D RngData;
uniform sampler2D RgbData;

varying vec2 vTexCoord;
uniform float SceneScale;


/////////////////////////////////////////////////////////////////////////
// Fresnel formulae
/////////////////////////////////////////////////////////////////////////

// N is outward normal (from medium to vacuum)
// Returns radiance gain on reflection (i.e., 1)
float Reflect(inout vec3 X, inout vec3 D, vec3 N)
{
	float cosi = dot(-D, N);
	bool entering = cosi > 0.0;
	if (!entering)
	{
		N *= -1.0; // Flip normal (so normal is always opposite to incident light direction)
	}
	// Reflect direction about normal, and displace ray start into reflected halfspace:
	float normalEpsilon = 2.0e-5*SceneScale;
	X += normalEpsilon*N;
	D -= 2.0*N*dot(D,N);
	return 1.0;
}


// D is direction of incident light
// N is normal pointing from medium to vacuum
float reflectionMetal(vec3 D, vec3 N, float ior, float k)
{
	float cosi = dot(-D, N);
	float cosip = abs(cosi);
	float cosi2 = cosip * cosip;
	float tmp = (ior*ior + k*k) * cosi2;
	float twoEtaCosi = 2.0 * ior * cosip;
	float Rparl2 = (tmp - twoEtaCosi + 1.0) / (tmp + twoEtaCosi + 1.0);
	float tmp_f = ior*ior + k*k;
	float Rperp2 = (tmp_f - twoEtaCosi + cosi2) / (tmp_f + twoEtaCosi + cosi2);
	return 0.5 * (Rparl2 + Rperp2);
}


// D is direction of incident light
// N is normal pointing from medium to vacuum
// returns radiance multiplier (varying for reflection, 0 for transmission, i.e. absorption)
float sampleMetal(inout vec3 X, inout vec3 D, vec3 N, float ior, float k, inout vec4 rnd)
{
	// @todo:  make sure that if light is emitted in interior of a metal, it terminates immediately.
	float R = reflectionMetal(D, N, ior, k);
	if (R >= rand(rnd)) // make reflectProb = R
	{
		// we must multiply subsequent radiance by factor (1.0 / reflectProb) to get correct counting
		// but as we chose reflectProb = R, this cancels with R, so that on reflection we should leave the radiance unchanged. 
		return Reflect(X, D, N);
	}
	else // absorption prob = (1-R)
	{
		return 0.0;
	}
}


// N is outward normal (from medium to vacuum).
// Returns radiance gain on transmission
float Transmit(inout vec3 X, inout vec3 D, vec3 N, float ior)
{
	// This applies of course only to dielectrics
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

	float r = ei/et;

	// Compute sint from Snells law:
	float sint = r * sqrt(max(0.0, 1.0 - cosi*cosi));

	// sint<=1.0 guaranteed as total internal reflection already handled
	float cost = sqrt(max(0.0, 1.0 - sint*sint)); 

	// Displace ray start into transmitted halfspace
	float normalEpsilon = 2.0e-5*SceneScale;
	X -= normalEpsilon*N;

	// Set transmitted direction
	D = r*D + (r*cosi - cost)*N; 

	// Transmitted radiance gets scaled by the square of the ratio of transmitted to incident IOR:
	return 1.0 / (r*r);
}


// D is direction of incident light
// N is normal pointing from medium to vacuum
float reflectionDielectric(vec3 D, vec3 N, float ior)
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
	float rParallel      = ( et*cosip - ei*cost) / ( et*cosip + ei*cost);
	float rPerpendicular = ( ei*cosip - et*cost) / ( ei*cosip + et*cost);
	return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);
}


// D is direction of incident light
// N is normal pointing from medium to vacuum
// returns radiance multiplier (1 for reflection, varying for transmission)
float sampleDielectric(inout vec3 X, inout vec3 D, vec3 N, float ior, inout vec4 rnd)
{
	float R = reflectionDielectric(D, N, ior);
	if (R >= rand(rnd)) // make reflectProb = R
	{
		// we must multiply subsequent radiance by factor (1.0 / reflectProb) to get correct counting
		// but as we chose reflectProb = R, this cancels with R, so that on reflection we should leave the radiance unchanged. 
		return Reflect(X, D, N);
	}
	else // refraction, prob = (1-R)
	{
		// we must multiply subsequent radiance by factor (1.0 / transmitProb) to get correct counting
		// but as we chose transmitProb = 1 - reflectProb = 1 - R, this cancels with (1-R) in the numerator, so that
		// on transmission, we should just multiply the radiance by the (et/ei)^2 gain factor (done inside transmission function)
		return Transmit(X, D, N, ior);
	}
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
	    SDF(X+eps.xyy) - SDF(X-eps.xyy),
	    SDF(X+eps.yxy) - SDF(X-eps.yxy),
	    SDF(X+eps.yyx) - SDF(X-eps.yyx) );
	return normalize(nor);
}

void raytrace(inout vec4 rnd, 
			  inout vec3 X, inout vec3 D,
			  inout vec3 rgb, float wavelength)
{
	if (length(rgb) < 1.0e-6) return;
	float minMarchDist = 1.0e-5*SceneScale;
	float maxMarchDist = 1.0e3*SceneScale;
	float t = 0.0;
	float h = 1.0;
    for( int i=0; i<MAX_MARCH_STEPS; i++ )
    {
		if (h<minMarchDist || t>maxMarchDist) break;
		h = abs(SDF(X + D*t));
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

	float wavelength = 360.0 + (750.0 - 360.0)*rgbw.w;
	raytrace(rnd, X, D, rgbw.rgb, wavelength);

	float sgn = sign( SDF(X) );

	gl_FragData[0] = vec4(X, sgn);
	gl_FragData[1] = vec4(D, 1.0);
	gl_FragData[2] = rnd;
	gl_FragData[3] = rgbw;
}


