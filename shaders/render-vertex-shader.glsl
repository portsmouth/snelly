
/////////////////////////////////////////////////
// Viewport vertex shader
/////////////////////////////////////////////////

// Line segment vertex positions come from a texture!
uniform sampler2D u_X;
attribute vec3 a_texCoord;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

void main()
{
	vec4 X = texture2D(u_X, a_texCoord.xy);
	vec3 P = mix(X.xyz, vec3(0.0), a_texCoord.z);

	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(P, 1.0);
}

