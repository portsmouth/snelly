
uniform sampler2D Radiance;
uniform sampler2D RngData;

// @todo:  camera details
uniform vec2 resolution;

uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;

uniform float camNear;
uniform float camFar;
const float camFovy; // degrees 


//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

SDF_FUNC

//////////////////////////////////////////////////////////////


bool hit(inout vec3 X, vec3 D)
{
	normalize(D);
	float minMarchDist = 1.0e-5*SceneScale;
	for (int i=0; i<MAX_MARCH_STEPS; i++)
	{
		float dist = abs(SDF(X));
		X += dist*D;
		if (dist < minMarchDist)
		{
			return true;
		}
		if (dist > 100.0*SceneScale)
		{
			return false;;
		}
	}
	return false;
}


vec3 localToWorld(vec3 N, vec3 wiL)
{
	vec3 T;
	if (abs(N.z) < abs(N.x))
	{
		T.x =  N.z;
		T.y =  0.0f;
		T.z = -N.x;
	}
	else
	{
		T.x =  0.0f;
		T.y =  N.z;
		T.z = -N.y;
	}
	vec3 B = N.cross(T);
	return T*wiL.x + B*wiL.y + N*wiL.z;
}


vec3 cosineSampleHemisphere(vec3 N, vec4 rnd)
{
	// sample disk
	float r = sqrtf(rand(rnd));
	float theta = 2.0 * M_PI * rand(rnd);
	vec2 p = vec2(r*cosf(theta), r*sinf(theta));

	// project
	float z = sqrtf(max(0.f, 1.f - p.x*p.x - p.y*p.y));
	return vec3(p.x, p.y, z);	
}


float computeClipDepth(float z, float zNear, float zFar)
{
	float zp = (zFar + zNear - 2.0f * zFar * zNear / z) / (zFar - zNear);
	zp = zp * 0.5f + 0.5f;
	return zp; // in [0,1] range as z ranges over [zNear, zFar]
}


void main()
{
	// Initialize world ray position
	vec3 X = camPos;

	// Compute world ray direction for this fragment
	vec2 ndc = -1.0 + 2.0* (gl_FragCoord.xy/resolution.xy);
	float aspect = resolution.x / resolution.y;
	float fh = 2.0*camNear*tan(0.5*radians(camFovy)); // frustum height
	float fw = aspect*fh;
	vec3 s = 0.5*(fw*ndc.x*camX + fh*ndc.y*camY);
	vec3 D = normalize(near*camDir + s); // ray direction

	// Raycast to first hit point
	vec4 rnd  = texture2D(RngData, vTexCoord);
	float zHit = camFar;
	vec3 L = vec3(0.0, 0.0, 0.0);
	
	if ( hit(X, D) )
	{
		zHit = length(X - camPos);

		// Construct a uniformly sampled AO 'shadow' ray in hemisphere of hit point
		vec3 N = NORMAL(X);

		// @todo: remind me why cosine-weighted sampling is the right thing here
		vec3 wiL = cosineSampleHemisphere(N, rnd);
		vec3 shadowRay = localToWorld(N, wiL);
		if ( !hit(X, shadowRay) )
		{
			// assume white background on ray escape
			L = vec3(1.0, 1.0, 1.0);
		}
	}

	// Write updated radiance and sample count
	vec4 oldL = texture2D(Radiance, vTexCoord);
	float oldN = oldL.w;
	float newN = oldN + 1.0;
	vec3 newL = (oldN*oldL + L) / newN;

	gl_FragColor = vec4(newL, newN);
	//gl_FragDepth = computeClipDepth(zHit, camNear, camFar);
}




