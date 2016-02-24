 
/////////////////////////////////////////////////
// Line vertex shader
/////////////////////////////////////////////////

precision highp float;

uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D RgbData;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

attribute vec3 TexCoord;
varying vec3 vColor;

void main()
{
	// Textures A and B contain line segment start and end points respectively
	vec3 posA = texture2D(PosDataA, TexCoord.xy).xyz;
	vec3 posB = texture2D(PosDataB, TexCoord.xy).xyz;

	// Line segment vertex position
	vec3 pos = mix(posA, posB, TexCoord.z);

	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(pos, 1.0);
	vColor = texture2D(RgbData, TexCoord.xy).rgb;
}

