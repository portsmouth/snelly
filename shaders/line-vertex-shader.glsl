
uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D RgbData;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;
uniform float sgn; // either +1 meaning we want to render only the exterior rays, 
                   //     or -1 meaning we want to render only the interior rays
                   //     or  0 meaning render all

attribute vec3 TexCoord;

varying vec3 vColor;

void main()
{
	// Textures A and B contain line segment start and end points respectively
	// (i.e. the geometry defined by this vertex shader is stored in textures)
	vec4 posA = texture2D(PosDataA, TexCoord.xy);
	vec4 posB = texture2D(PosDataB, TexCoord.xy);

	float sgnA = posA.w; // SDF sign: +1.0 for exterior rays, or -1.0 for interior rays
	float kill = (1.0-abs(sgn)) + abs(sgn)*0.5*abs(sgnA + sgn); 
	
	// Line segment vertex position (either posA or posB)
	vec3 pos = mix(posA.xyz, posB.xyz, TexCoord.z);

	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(pos, 1.0);
	vColor = kill * texture2D(RgbData, TexCoord.xy).rgb;
}

