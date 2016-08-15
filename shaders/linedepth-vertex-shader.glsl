
uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D Depth;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

attribute vec3 TexCoord;

//varying float zCoord;
varying vec2 vTexCoord;

varying float eye_z;
varying float eye_w;


void main()
{
	// Textures A and B contain line segment start and end points respectively
	// (i.e. the geometry defined by this vertex shader is stored in textures)
	vec3 posA = texture2D(PosDataA, TexCoord.xy).xyz;
	vec3 posB = texture2D(PosDataB, TexCoord.xy).xyz;

	// Line segment vertex position (either posA or posB)
	vec3 pos = mix(posA, posB, TexCoord.z);

	vec4 pEye = u_modelViewMatrix * vec4(pos, 1.0);
	eye_z = -pEye.z; // Keep eye space depth for fragment shader

	gl_Position = u_projectionMatrix * pEye; // Transform to clip space for vertex shader
}

