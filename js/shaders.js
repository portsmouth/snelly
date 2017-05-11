var Shaders = {

'comp-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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
uniform bool drawExterior;

varying vec2 vTexCoord;

void main() 
{
	vec3 fluence = vec3(0.0, 0.0, 0.0);
	if (drawInterior)
	{
		fluence += texture2D(FluenceInt, vTexCoord).rgb;
	}
	if (drawExterior)
	{
		fluence += texture2D(FluenceExt, vTexCoord).rgb;
	}

	vec3 phi = invNumPaths * fluence; // normalized fluence

	// Apply exposure 
	float gain = pow(2.0, exposure);
	float r = gain*phi.x; 
	float g = gain*phi.y; 
	float b = gain*phi.z;
	
	// Reinhard tonemap
	vec3 C = vec3(r/(1.0+r), g/(1.0+g), b/(1.0+b)) ;

	// Apply gamma
	C = pow(C, vec3(invGamma));

	gl_FragColor = vec4(C, 1.0);
}
`,

'comp-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

'filter-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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
uniform vec2 resolution;

#define RADIANCE_TOLERANCE 1.0e-6


float filter(float r2)
{
    return max(0.0, 1.0-r2);
}

void main()
{
    vec2 pixel = gl_FragCoord.xy;

    // Average current radiance over sample set within filter radius
    float maxx = resolution.x-1.0;
    float maxy = resolution.y-1.0;
    vec2 invRes = 1.0/resolution.xy;
    float invWidth = 0.5/float(DOWN_RES);

    vec3 mean = vec3(0.0);
    float norm = 0.0;
    for (int i=-2*DOWN_RES; i<=2*DOWN_RES; ++i)
    {
        float _i = max(min(maxx, pixel.x+float(i)), 0.0);
        float u = _i*invRes.x;
        float dx = float(i)*invWidth;
        float dx2 = dx*dx;

        for (int j=-2*DOWN_RES; j<2*DOWN_RES; ++j)
        {
            float _j = max(min(maxy, pixel.y+float(j)), 0.0);
            float v = _j*invRes.y;
            float dy = float(j)*invWidth;
            
            float r2 = dx2 + dy*dy;
            float f = filter(r2);
            if (f>1.0e-3)
            {
                vec4 L = texture2D(Radiance, vec2(u,v));  
                float weight = f * L.w;
                mean += weight * L.xyz;
                norm += weight;
            }
        }
    }

    mean /= max(norm, RADIANCE_TOLERANCE);
    gl_FragData[0] = vec4(mean, 1.0);
}
`,

'filter-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

'init-fragment-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

	// Sample photon wavelength via the inverse CDF of the emission spectrum
	// (here w is spectral offset, i.e. wavelength = 360.0 + (750.0 - 360.0)*w)
	// (linear interpolation into the inverse CDF texture and RGB table should ensure smooth sampling over the range)
    float w = texture2D(ICDF, vec2(rand(seed), 0.5)).r;
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
	u = normalize(u);
	vec3 v = cross(EmitterDir, u);
	vec3 pos = EmitterPos + rPos*(u*cos(phiPos) + v*sin(phiPos)); 

	// Emit in a cone with the given spread
	float spreadAngle = 0.5*abs(EmitterSpread)*M_PI/180.0;
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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


vec3 xyzToRgb(vec3 XYZ)
{
    // (using sRGB color space)
    vec3 RGB;
    RGB.r =  3.2404542*XYZ.x - 1.5371385*XYZ.y - 0.4985314*XYZ.z;
    RGB.g = -0.9692660*XYZ.x + 1.8760108*XYZ.y + 0.0415560*XYZ.z;
    RGB.b =  0.0556434*XYZ.x - 0.2040259*XYZ.y + 1.0572252*XYZ.z;
    return RGB;
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
                vec3 RGB = xyzToRgb(XYZ);
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
`,

'pathtracer-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/*
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
*/

void main()
{
	// @todo ...

	/*
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
	*/


	float dist = -1.0;
	gl_FragColor = encode_float(dist);
}
`,

'pick-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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
varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float whitepoint;

void main()
{
	vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
	float X = L.x;
	float Y = L.y;
	float Z = L.z;
	float sum = X + Y + Z;
	float x = X / sum;
	float y = Y / sum; 

	// compute Reinhard tonemapping scale factor
	float scale = (1.0 + Y/(whitepoint*whitepoint)) / (1.0 + Y);
	Y *= scale;
	X = x * Y / y;
	Z = (1.0 - x - y) * (Y / y);

	// Convert XYZ tristimulus to sRGB color space
	float r =  3.2406*X - 1.5372*Y - 0.4986*Z;
	float g = -0.9689*X + 1.8758*Y + 0.0415*Z;
	float b =  0.0557*X - 0.2040*Y + 1.0570*Z;

	vec3 Lp = vec3(r, g, b);
	vec3 S = pow(Lp, vec3(invGamma));

	gl_FragColor =vec4(S, 0.0);
}
`,

'tonemapper-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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
`,

'trace-vertex-shader': `
#extension GL_EXT_draw_buffers : require
//#extension GL_EXT_frag_depth : require
precision highp float;

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
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


float saturate(float x) { return max(0.0, min(1.0, x)); }

// clip space depth calculation, given eye space depth
float computeClipDepth(float ze, float zNear, float zFar)
{
	float zn = (zFar + zNear - 2.0*zFar*zNear/ze) / (zFar - zNear); // compute NDC depth
	float zb = zn*0.5 + 0.5;                                        // convert to clip depth
	return saturate(zb); // in [0,1] range as z ranges over [zNear, zFar]
}


///
/// A neat trick to return a float value from a webGL fragment shader
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