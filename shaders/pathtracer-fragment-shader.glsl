
uniform sampler2D Radiance;
uniform sampler2D RngData;
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


//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

SDF_FUNC

//////////////////////////////////////////////////////////////


bool hit(inout vec3 X, vec3 D, inout int numSteps)
{
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
	if (t<maxMarchDist) return true;
	return false;
}


vec3 NORMAL( in vec3 pos )
{
	// Compute normal as gradient of SDF
	float normalEpsilon = 2.0e-5*SceneScale;
	vec3 eps = vec3(normalEpsilon, 0.0, 0.0);
	vec3 nor = vec3(
	    SDF(pos+eps.xyy) - SDF(pos-eps.xyy),
	    SDF(pos+eps.yxy) - SDF(pos-eps.yxy),
	    SDF(pos+eps.yyx) - SDF(pos-eps.yyx) );
	return normalize(nor);
}


vec3 localToWorld(vec3 N, vec3 wiL)
{
	vec3 T;
	if (abs(N.z) < abs(N.x))
	{
		T.x =  N.z;
		T.y =  0.0;
		T.z = -N.x;
	}
	else
	{
		T.x =  0.0;
		T.y =  N.z;
		T.z = -N.y;
	}
	vec3 B = cross(N, T);
	return T*wiL.x + B*wiL.y + N*wiL.z;
}


vec3 cosineSampleHemisphere(vec3 N, inout vec4 rnd)
{
	// sample disk
	float r = sqrt(rand(rnd));
	float theta = 2.0*M_PI*rand(rnd);
	vec2 p = vec2(r*cos(theta), r*sin(theta));

	// project
	float z = sqrt(max(0.0, 1.0 - p.x*p.x - p.y*p.y));
	return vec3(p.x, p.y, z);	
}

/*
float computeClipDepth(float z, float zNear, float zFar)
{
	float zp = (zFar + zNear - 2.0*zFar*zNear/z) / (zFar - zNear);
	zp = zp * 0.5 + 0.5;
	return zp; // in [0,1] range as z ranges over [zNear, zFar]
}
*/

void main()
{
	vec4 rnd = texture2D(RngData, vTexCoord);

	// Initialize world ray position
	vec3 X = camPos;

	// Jitter over pixel
	vec2 pixel = gl_FragCoord.xy;
	pixel += -0.5 + 0.5*vec2(rand(rnd), rand(rnd));

	// Compute world ray direction for this fragment
	vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
	float fh = camNear*tan(0.5*radians(camFovy)) / camZoom; // frustum height
	float fw = camAspect*fh;
	vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
	vec3 D = normalize(camNear*camDir + s); // ray direction

	// Raycast to first hit point
	//float zHit = camFar;
	vec3 L = vec3(0.0, 0.0, 0.0);
	int numSteps;	
	if ( hit(X, D, numSteps) )
	{
		//zHit = length(X - camPos);
		vec3 N = NORMAL(X);
		L = 0.5*(1.0+N);
	}

	// Write updated radiance and sample count
	vec4 oldL = texture2D(Radiance, vTexCoord);
	float oldN = oldL.w;
	float newN = oldN + 1.0;
	vec3 newL = (oldN*oldL.rgb + L) / newN;

	gl_FragData[0] = vec4(newL, newN); 
	gl_FragData[1] = rnd;

	//gl_FragDepth = computeClipDepth(zHit, camNear, camFar);
}




