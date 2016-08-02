var Shaders = {

'comp-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

uniform sampler2D Frame;
uniform float Exposure;
varying vec2 vTexCoord;

void main() 
{
	// @todo: expose gamma here in UI.
	// @todo: 'proper' tonemapping here.
	gl_FragColor = vec4(pow(texture2D(Frame, vTexCoord).rgb*Exposure, vec3(1.0/2.2)), 1.0);
}
`,

'comp-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

attribute vec3 Position;
attribute vec2 TexCoord;

varying vec2 vTexCoord;

void main(void)
{
	gl_Position = vec4(Position, 1.0);
	vTexCoord = TexCoord;
}
`,

'init-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

uniform sampler2D RngData;
uniform sampler2D WavelengthToRgb;
uniform sampler2D ICDF;

uniform vec3 EmitterPos;
uniform vec3 EmitterDir;
uniform float EmitterRadius;
uniform float EmitterSpread; // in degrees
uniform float EmitterPower;

varying vec2 vTexCoord;

void main()
{
	vec4 seed = texture2D(RngData, vTexCoord);

	// Sample photon wavelength from CDF of emission spectrum
	// (here w is spectral offset, i.e. wavelength = 360.0 + (750.0 - 360.0)*w)
    float w = texture2D(ICDF, vec2(rand(seed), 0.5)).r + rand(seed)*(1.0/256.0);
  	vec3 rgb = EmitterPower * texture2D(WavelengthToRgb, vec2(w, 0.5)).rgb;

	// Make emission cross-section circular
	float rPos   = EmitterRadius*sqrt(rand(seed));
	float phiPos = 2.0*M_PI*rand(seed);
	vec3 X = vec3(1.0, 0.0, 0.0);
	vec3 Z = vec3(0.0, 0.0, 1.0);
	vec3 u = cross(Z, EmitterDir);
	if ( length(u) < 1.0e-3 )
	{
		u = cross(X, EmitterDir);
	}
	vec3 v = cross(EmitterDir, u);
	vec3 pos = EmitterPos + rPos*(u*cos(phiPos) + v*sin(phiPos)); 

	// Emit in a cone with the given spread
	float spreadAngle = abs(EmitterSpread)*M_PI/180.0;
	float rDir = min(tan(spreadAngle), 1.0e6) * sqrt(rand(seed));
	float phiDir = 2.0*M_PI*rand(seed);
	vec3 dir = normalize(EmitterDir + rDir*(u*cos(phiDir) + v*sin(phiDir)));
	
	gl_FragData[0] = vec4(pos, 1.0);
	gl_FragData[1] = vec4(dir, 1.0);
	gl_FragData[2] = seed;
	gl_FragData[3] = vec4(rgb, w);
}
`,

'init-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

attribute vec3 Position;
attribute vec2 TexCoord;

varying vec2 vTexCoord;

void main() 
{
	gl_Position = vec4(Position, 1.0);
	vTexCoord = TexCoord;
}
`,

'line-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

varying vec3 vColor;

void main() 
{
	gl_FragColor = vec4(vColor, 1.0);
}
`,

'line-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D RgbData;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

attribute vec3 TexCoord;
varying vec3 vColor;

void main()
{
	// Textures A and B contain line segment start and end points respectively
	// (i.e. the geometry defined by this vertex shader comes from textures!)
	vec3 posA = texture2D(PosDataA, TexCoord.xy).xyz;
	vec3 posB = texture2D(PosDataB, TexCoord.xy).xyz;

	// Line segment vertex position
	vec3 pos = mix(posA, posB, TexCoord.z);

	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(pos, 1.0);
	vColor = texture2D(RgbData, TexCoord.xy).rgb;
}
`,

'pass-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

uniform sampler2D Frame;
varying vec2 vTexCoord;

void main() 
{
	gl_FragColor = vec4(texture2D(Frame, vTexCoord).rgb, 1.0);
}
`,

'pass-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

attribute vec3 Position;
attribute vec2 TexCoord;

varying vec2 vTexCoord;

void main(void)
{
	gl_Position = vec4(Position, 1.0);
	vTexCoord = TexCoord;
}
`,

'pathtracer-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

uniform sampler2D Radiance;
uniform sampler2D RngData;
varying vec2 vTexCoord;

// @todo:  camera details
uniform vec2 resolution;

uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;

uniform float camNear;
uniform float camFar;
uniform float camFovy; // degrees 

uniform float SceneScale;


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


vec3 cosineSampleHemisphere(vec3 N, vec4 rnd)
{
	// sample disk
	float r = sqrt(rand(rnd));
	float theta = 2.0 * M_PI * rand(rnd);
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
	// Initialize world ray position
	vec3 X = camPos;

	// Compute world ray direction for this fragment
	vec2 ndc = -1.0 + 2.0* (gl_FragCoord.xy/resolution.xy);

	/*
	float aspect = resolution.x / resolution.y;
	float fh = 2.0*camNear*tan(0.5*radians(camFovy)); // frustum height
	float fw = aspect*fh;
	vec3 s = 0.5*(fw*ndc.x*camX + fh*ndc.y*camY);
	vec3 D = normalize(camNear*camDir + s); // ray direction

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
	vec3 newL = (oldN*oldL.rgb + L) / newN;
	*/

	//vec4 oldL = texture2D(Radiance, vTexCoord);

	vec3 L = vec3(ndc.x, ndc.y, 0.0);
	gl_FragData[0] = vec4(L, 1.0);
	
	vec4 rnd  = texture2D(RngData, vTexCoord);
	rand(rnd);
	gl_FragData[1] = rnd;

	//gl_FragDepth = computeClipDepth(zHit, camNear, camFar);
}
`,

'pathtracer-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

uniform sampler2D Radiance;
varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;

void main()
{
	vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
	float r = L.x; 
	float g = L.y; 
	float b = L.z;
	vec3 Lp = vec3(r/(1.0+r), g/(1.0+g), b/(1.0+b));

	gl_FragColor = vec4(pow(Lp, vec3(invGamma)), 1.0);
}
`,

'tonemapper-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

attribute vec3 Position;
attribute vec2 TexCoord;
varying vec2 vTexCoord;

void main() 
{
	gl_Position = vec4(Position, 1.0);
	vTexCoord = TexCoord;
}
`,

'trace-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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

void raytrace(inout vec4 rnd, 
			  inout vec3 X, inout vec3 D,
			  inout vec3 rgb, float wavelength)
{
	if (length(rgb) < 1.0e-6) return;

	bool hit = false;
	normalize(D);
	float minMarchDist = 1.0e-5*SceneScale;

	for (int i=0; i<MAX_MARCH_STEPS; i++)
	{
		float dist = abs(SDF(X));
		X += dist*D;
		if (dist < minMarchDist)
		{
			hit = true;
			break;
		}
		if (dist > 1000.0*SceneScale)
		{
			break;
		}
	}

	if (!hit)
	{
		X += 1000.0*SceneScale*D;
		rgb *= 0.0; // terminate ray
	}
	else
	{
		rgb *= SAMPLE(X, D, NORMAL(X), wavelength, rnd);
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

	gl_FragData[0] = vec4(X, 1.0);
	gl_FragData[1] = vec4(D, 1.0);
	gl_FragData[2] = rnd;
	gl_FragData[3] = rgbw;
}
`,

'trace-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL float point pseudorandom number generator, from
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