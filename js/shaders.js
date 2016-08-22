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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D FluenceInt;
uniform sampler2D FluenceExt;

uniform float invNumPaths;
uniform float exposure;
uniform float invGamma;
uniform bool drawInterior;

varying vec2 vTexCoord;

void main() 
{
	vec3 fluence = texture2D(FluenceExt, vTexCoord).rgb;
	if (drawInterior)
	{
		fluence += texture2D(FluenceInt, vTexCoord).rgb;
	}

	vec3 L = invNumPaths * pow(10.0, exposure) * fluence;
	float r = L.x; 
	float g = L.y; 
	float b = L.z;
	vec3 Lp = vec3(r/(1.0+r), g/(1.0+g), b/(1.0+b));
	vec3 S = pow(Lp, vec3(invGamma));

	gl_FragColor = vec4(S, 1.0);
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D RngData;
uniform sampler2D WavelengthToRgb;
uniform sampler2D ICDF;

uniform vec3 EmitterPos;
uniform vec3 EmitterDir;
uniform float EmitterRadius;
uniform float EmitterSpread; // in degrees

varying vec2 vTexCoord;

void main()
{
	vec4 seed = texture2D(RngData, vTexCoord);

	// Sample photon wavelength from CDF of emission spectrum
	// (here w is spectral offset, i.e. wavelength = 360.0 + (750.0 - 360.0)*w)
    float w = texture2D(ICDF, vec2(rand(seed), 0.5)).r;// + rand(seed)*(1.0/256.0);
  	vec3 rgb = texture2D(WavelengthToRgb, vec2(w, 0.5)).rgb;

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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D RgbData;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;
uniform float sgn; // either +1 meaning we want to render only the exterior rays, 
                   //     or -1 meaning we want to render only the interior rays
                   //     or  0 meaning render all

attribute vec3 TexCoord;

varying vec3 vColor;

void main()
{
	// Textures A and B contain line segment start and end points respectively
	// (i.e. the geometry defined by this vertex shader is stored in textures)
	vec4 posA = texture2D(PosDataA, TexCoord.xy);
	vec4 posB = texture2D(PosDataB, TexCoord.xy);

	float sgnA = posA.w; // SDF sign: +1.0 for exterior rays, or -1.0 for interior rays
	float kill = (1.0-abs(sgn)) + abs(sgn)*0.5*abs(sgnA + sgn); 
	
	// Line segment vertex position (either posA or posB)
	vec3 pos = mix(posA.xyz, posB.xyz, TexCoord.z);

	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(pos, 1.0);
	vColor = kill * texture2D(RgbData, TexCoord.xy).rgb;
}
`,

'linedepth-fragment-shader': `
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D Depth;

uniform float camNear;
uniform float camFar;

varying float eye_z;
uniform vec2 resolution;


void main() 
{
	vec2 viewportTexCoords = gl_FragCoord.xy/resolution;
	float destClipDepth = unpack_depth( texture2D(Depth, viewportTexCoords) );
	
	float sourceClipDepth = computeClipDepth(eye_z, camNear, camFar);
	if (sourceClipDepth < destClipDepth)
	{
		gl_FragColor = pack_depth(sourceClipDepth);
	}
	else
	{
		gl_FragColor = pack_depth(destClipDepth);
	}
}
`,

'linedepth-vertex-shader': `
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D Depth;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

attribute vec3 TexCoord;

varying float eye_z;

void main()
{
	// Textures A and B contain line segment start and end points respectively
	// (i.e. the geometry defined by this vertex shader is stored in textures)
	vec3 posA = texture2D(PosDataA, TexCoord.xy).xyz;
	vec3 posB = texture2D(PosDataB, TexCoord.xy).xyz;

	// Line segment vertex position (either posA or posB)
	vec3 pos = mix(posA, posB, TexCoord.z);

	vec4 pEye = u_modelViewMatrix * vec4(pos, 1.0);
	eye_z = -pEye.z; // Keep eye space depth for fragment shader

	gl_Position = u_projectionMatrix * pEye; // Transform to clip space for vertex shader
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D WaveBuffer;
varying vec2 vTexCoord;

void main() 
{
	gl_FragColor = vec4(texture2D(WaveBuffer, vTexCoord).rgba);
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D Radiance;
uniform sampler2D RngData;
uniform sampler2D Depth;

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

LIGHTING_FUNC

//////////////////////////////////////////////////////////////


bool hit(inout vec3 X, vec3 D, inout int numSteps)
{
	float minMarchDist = 1.0e-5*SceneScale;
	float maxMarchDist = 1.0e5*SceneScale;
	float t = 0.0;
	float h = 1.0;
    for( int i=0; i<MAX_MARCH_STEPS; i++ )
    {
		if (h<minMarchDist || t>maxMarchDist) break;
		h = SDF(X + D*t);
        t += h;
    }
    X += t*D;
	if (t<maxMarchDist) return true;
	return false;
}


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
	float zEye = camFar;
	vec3 L = vec3(0.0, 0.0, 0.0);
	int numSteps;	
	if ( hit(X, D, numSteps) )
	{
		zEye = dot(X - camPos, camDir);
		vec3 N = NORMAL(X);
		vec3 V = normalize(camPos-X);
		L = LIGHTING(V, N);
	}

	float clipDepth = computeClipDepth(zEye, camNear, camFar);

	// Write updated radiance and sample count
	vec4 oldL = texture2D(Radiance, vTexCoord);
	float oldN = oldL.w;
	float newN = oldN + 1.0;
	vec3 newL = (oldN*oldL.rgb + L) / newN;

	gl_FragData[0] = vec4(newL, newN);
	gl_FragData[1] = rnd;
	gl_FragData[2] = pack_depth(clipDepth);
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

'pick-fragment-shader': `
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform float ndcX; // NDC coordinates of pick
uniform float ndcY;

uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;

uniform float camNear;
uniform float camFovy; // degrees 
uniform float camZoom;
uniform float camAspect;

uniform float SceneScale;


//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

SDF_FUNC

//////////////////////////////////////////////////////////////


bool hit(inout vec3 X, vec3 D)
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

void main()
{
	// Initialize world ray position
	vec3 X = camPos;

	// Compute world ray direction for this fragment
	vec2 ndc = vec2(ndcX, ndcY);
	float fh = camNear*tan(0.5*radians(camFovy)) / camZoom;
	float fw = camAspect*fh;
	vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
	vec3 D = normalize(camNear*camDir + s); // ray direction

	// Raycast to first hit point
	float dist;
	if ( hit(X, D) )
	{	
		dist = length(X - camPos);
	}
	else
	{
		dist = -1.0;
	}
	
	gl_FragColor = encode_float(dist);
}
`,

'pick-vertex-shader': `
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
}

uniform sampler2D Radiance;
uniform sampler2D DepthSurface;
uniform sampler2D DepthLight;

varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float alpha;
uniform bool enableDepthTest;


void main()
{
	vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
	float r = L.x; 
	float g = L.y; 
	float b = L.z;
	vec3 Lp = vec3(r/(1.0+r), g/(1.0+g), b/(1.0+b));
	vec3 S = pow(Lp, vec3(invGamma));
	vec3 Sp = S * alpha;

	float surfaceDepth = unpack_depth(texture2D(DepthSurface, vTexCoord));
	float   lightDepth = unpack_depth(texture2D(DepthLight, vTexCoord));

	// Composite surface with light ray fragment
	float A = 0.0;
	if (enableDepthTest && (lightDepth > surfaceDepth))
	{
		// light behind surface
		A = alpha;
	}
	
	gl_FragColor = vec4(Sp, A);
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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
	float maxMarchDist = 1.0e5*SceneScale;
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

/// Distance field utilities

// Union
float opU( float d1, float d2 ) { return min(d1,d2); }

// Subtraction
float opS(float A, float B) { return max(-B, A); }

// Intersection
float opI( float d1, float d2 ) { return max(d1,d2); }


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a extract a float value from a fragment shader
///
float shift_right (float v, float amt) { 
    v = floor(v) + 0.5; 
    return floor(v / exp2(amt)); 
}
float shift_left (float v, float amt) { 
    return floor(v * exp2(amt) + 0.5); 
}
float mask_last (float v, float bits) { 
    return mod(v, shift_left(1.0, bits)); 
}
float extract_bits (float num, float from, float to) { 
    from = floor(from + 0.5); to = floor(to + 0.5); 
    return mask_last(shift_right(num, from), to - from); 
}
vec4 encode_float (float val) { 
    if (val == 0.0) return vec4(0, 0, 0, 0); 
    float sign = val > 0.0 ? 0.0 : 1.0; 
    val = abs(val); 
    float exponent = floor(log2(val)); 
    float biased_exponent = exponent + 127.0; 
    float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0; 
    float t = biased_exponent / 2.0; 
    float last_bit_of_biased_exponent = fract(t) * 2.0; 
    float remaining_bits_of_biased_exponent = floor(t); 
    float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0; 
    float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0; 
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0; 
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0; 
    return vec4(byte4, byte3, byte2, byte1); 
}


vec4 pack_depth(const in float depth)
{
    return vec4(depth, 0.0, 0.0, 1.0);
}

float unpack_depth(const in vec4 rgba_depth)
{
    return rgba_depth.r;
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