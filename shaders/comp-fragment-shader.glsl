
uniform sampler2D Fluence;
uniform float invNumPaths;
uniform float exposure;
uniform float invGamma;
varying vec2 vTexCoord;

void main() 
{
	vec3 L = invNumPaths * pow(10.0, exposure) * texture2D(Fluence, vTexCoord).rgb;
	gl_FragColor = vec4(pow(L, vec3(invGamma)), 1.0);
}
