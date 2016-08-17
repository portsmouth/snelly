
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