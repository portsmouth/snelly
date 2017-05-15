
uniform sampler2D Radiance;
varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float whitepoint;


void constrain_rgb(inout vec3 RGB)
{
	float w;
	w = (0.0 < RGB.r) ? 0.0 : RGB.r;
	w = (w   < RGB.g) ? 0.0 : RGB.g;
	w = (w   < RGB.b) ? 0.0 : RGB.b;
	w = -w;
	if (w>0.0)
	{
		RGB.r += w; RGB.g += w; RGB.b += w;
	}
}

void main()
{
	vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
	//vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
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

	// convert XYZ tristimulus to sRGB color space
	vec3 RGB;
	RGB.r =  3.2406*X - 1.5372*Y - 0.4986*Z;
	RGB.g = -0.9689*X + 1.8758*Y + 0.0415*Z;
	RGB.b =  0.0557*X - 0.2040*Y + 1.0570*Z;

	// deal with out-of-gamut RGB.
	constrain_rgb(RGB);

	// apply gamma correction
	vec3 S = pow(abs(RGB), vec3(invGamma));

	gl_FragColor =vec4(S, 0.0);
}