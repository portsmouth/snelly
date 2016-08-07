
uniform sampler2D Radiance;
varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float alpha;

void main()
{
	vec3 L = exposure * texture2D(Radiance, vTexCoord).rgb;
	float r = L.x; 
	float g = L.y; 
	float b = L.z;
	vec3 Lp = vec3(r/(1.0+r), g/(1.0+g), b/(1.0+b));

	gl_FragColor = vec4(pow(Lp, vec3(invGamma)), alpha);
}