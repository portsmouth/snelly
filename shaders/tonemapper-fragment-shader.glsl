
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