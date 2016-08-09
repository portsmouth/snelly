
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

///
/// A neat trick to encode a float value in the frag color
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




