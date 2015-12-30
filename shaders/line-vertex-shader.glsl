
/////////////////////////////////////////////////
// Line vertex shader
/////////////////////////////////////////////////

// Line segment vertex positions come from a texture!
uniform sampler2D u_X;
attribute vec3 TexCoord;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

void main()
{
	vec4 X = texture2D(u_X, TexCoord.xy);
	vec3 P = mix(X.xyz, vec3(0.0), TexCoord.z);
	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(P, 1.0);
}



/* // debug version
attribute vec3 a_position;
attribute vec2 a_texCoord;
varying vec2 v_texCoord;
void main()
{
	gl_Position = vec4(a_position, 1.0);
	v_texCoord = a_texCoord;
}
*/
