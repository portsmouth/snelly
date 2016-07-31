var Shaders = {

'comp-fragment-shader': `
#extension GL_EXT_draw_buffers : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
	// @todo: expose gamma here in UI
	gl_FragColor = vec4(pow(texture2D(Frame, vTexCoord).rgb*Exposure, vec3(1.0/2.2)), 1.0);
}
`,

'comp-vertex-shader': `
#extension GL_EXT_draw_buffers : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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

'trace-fragment-shader': `
#extension GL_EXT_draw_buffers : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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
	bool hit = false;
	normalize(D);
	float minMarchDist = 1.0e-5*SceneScale;

	//for (int i=0; i<MAX_MARCH_STEPS; i++)
	for (int i=0; i<256; i++)
	{
		float dist = abs(SDF(X));
		X += dist*D;
		if (dist < minMarchDist)
		{
			hit = true;
			break;
		}
		if (dist > 100.0*SceneScale)
		{
			break;
		}
	}

	if (!hit)
	{
		X += SceneScale*D;
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
precision highp float;

#define M_PI 3.1415926535897932384626433832795


/// @todo: where does this come from? What is it?

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