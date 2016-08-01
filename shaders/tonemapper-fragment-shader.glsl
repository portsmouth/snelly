
uniform sampler2D samplerHDR;
uniform float exposure;
uniform float invGamma;
in vec2 varTexCoord0;

void main()
{
	vec3 L = exposure * texture(samplerHDR, varTexCoord0).rgb;
	
	float r = L.x; 
	float g = L.y; 
	float b = L.z;
	
	vec3 Lp = vec3(r/(1.0+r), g/(1.0+g), b/(1.0+b));

	outColor = vec4(pow(Lp, vec3(invGamma)), 1.0f);
}