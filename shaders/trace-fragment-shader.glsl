
uniform sampler2D PosData;
uniform sampler2D DirData;
uniform sampler2D RngData;
uniform sampler2D RgbData;

varying vec2 vTexCoord;

uniform float SceneScale;
uniform int MaxMarchSteps;


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
	float normalEpsilon = 5.0e-5 * SceneScale;
	X += normalEpsilon*N;
	D -= 2.0*N*dot(D,N);
	return 1.0;
}


//////////////////////////////////////////////////////////////
// Dielectric formulae
//////////////////////////////////////////////////////////////


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
	float normalEpsilon = 5.0e-5 * SceneScale;
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
float sampleDielectric(inout vec3 X, inout vec3 D, vec3 N, float ior)
{
	float R = reflectionDielectric(D, N, ior);
	float rnd = rand(state);
	if (R >= rnd) // make reflectProb = R
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
		return Transmit(X, D, N, ior, gain);
	}
}


//////////////////////////////////////////////////////////////
// Metal formulae
//////////////////////////////////////////////////////////////

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
float sampleMetal(inout vec3 X, inout vec3 D, vec3 N, float ior, float k)
{
	float R = reflectionMetal(D, N, ior, k);
	if (R >= rnd) // make reflectProb = R
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




//////////////////////////////////////////////////////////////
// @todo: paste this code in dynamically, based on current scene
//////////////////////////////////////////////////////////////

uniform float _radius;                
float SDF(vec3 X)                     
{                                     
	return length(X) - _radius;       
}                                     


//////////////////////////////////////////////////////////////
// @todo: paste this code in dynamically, based on current material
//////////////////////////////////////////////////////////////

float ior(float lnm) 
{
	return 1.5;
}

float sample(inout vec3 X, inout vec3 D, vec3 N, float lnm)
{
	return sampleDielectric(X, D, N, ior(lnm));
}



//////////////////////////////////////////////////////////////
// Main SDF tracing loop
//////////////////////////////////////////////////////////////

vec3 calcNormal(in vec3 X)
{
	// Compute normal as gradient of SDF
	float eps = 1.0e-5 * SceneScale;
	vec3 N = vec3( SDF(X+eps) - SDF(X-eps),
				   SDF(X+eps) - SDF(X-eps),
				   SDF(X+eps) - SDF(X-eps) );
	return normalize(N);
}


void raytrace(inout vec3 X, inout vec3 D,
			  inout vec4 rgbLambda, 
			  inout vec4 state)
{
	const float radianceEpsilon = 1.0e-7;
	if ( length(rgbLambda.rgb) < radianceEpsilon ) return;

	float totalDist = 0.0;
	bool hit = false;

	float minMarchDist = 1.0e-5 * SceneScale;
	for (int i=0; i<MaxMarchSteps; i++)
	{
		X += totalDist*D;
		float dist = abs( SDF(X) );
		totalDist += dist;
		if (dist < minMarchDist)
		{
			hit = true;
			break;
		}
	}

	if (!hit)
	{
		X += maxDist*D;
		rgbLambda.rgb *= 0.0; // terminate ray
	}
	else
	{
		// Hit the surface. Calculate normal there:
		vec3 N = calcNormal(Xp);

		// Sample new direction, and update radiance accordingly:
		float lambda = rgbLambda.w;
		rgbLambda.rgb *= sample(X, D, N, lambda);
	}
}


void main()
{
	vec3 X         = texture2D(PosData, vTexCoord).xyz;
	vec3 D         = texture2D(DirData, vTexCoord).xyz;
	vec4 state     = texture2D(RngData, vTexCoord);
	vec4 rgbLambda = texture2D(RgbData, vTexCoord);

	raytrace(X, D, rgbLambda, state);

	gl_FragData[0] = vec4(X, 1.0);
	gl_FragData[1] = vec4(D, 1.0);
	gl_FragData[2] = state;
	gl_FragData[3] = rgbLambda;
}


