
uniform sampler2D Radiance;
uniform sampler2D DepthSurface;
uniform sampler2D DepthLight;

varying vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float alpha;
uniform int enableDepthTest;


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
	float A;
	if (lightDepth < surfaceDepth)
	{
		// light in front of surface
		A = 0.0;

		// with gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
		// this should produce
		//    rgb = Sp + Light
		//      a = A + 1.0*(1-A) = 1.0
    }
    else
    {
    	// light behind surface
    	A = alpha;

		// with gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
		// this should produce
		//    rgb = Sp + (1-alpha)*Light
		//      a = A + 1.0*(1-A) = 1.0
    }

	gl_FragColor = vec4(Sp, A);



	//gl_FragColor = vec4(packedLightClipDepth.rgb, 1);
	//gl_FragColor = vec4(lightClipDepth, lightClipDepth, lightClipDepth, 1);
}