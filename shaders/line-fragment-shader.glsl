
/////////////////////////////////////////////////
// Viewport fragment shader
/////////////////////////////////////////////////

precision mediump float;

uniform sampler2D u_X;
//varying vec2 v_texCoord;

void main() 
{
	//vec4 X = texture2D(u_X, v_texCoord);
	//gl_FragColor = vec4(X.xyz, 1.0);
	gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}
