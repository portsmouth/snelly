
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
