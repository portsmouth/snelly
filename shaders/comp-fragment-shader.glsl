
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
