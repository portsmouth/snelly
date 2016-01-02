
/////////////////////////////////////////////////
// Raytrace fragment shader
/////////////////////////////////////////////////

#extension GL_EXT_draw_buffers : require
precision mediump float;

uniform sampler2D PosData;
uniform sampler2D DirData;
uniform sampler2D RngData;
uniform sampler2D RgbData;

varying vec2 vTexCoord;

void main() 
{
	vec3 X         = texture2D(PosData, vTexCoord).xyz;
	vec3 D         = texture2D(DirData, vTexCoord).xyz;
	vec4 state     = texture2D(RngData, vTexCoord);
	vec4 rgbLambda = texture2D(RgbData, vTexCoord);

	// Read X, W, F
	// Do:  X' = trace(X, W) given F
	// Write X' into fragment (of ray buffer)
	vec3 Xp = X + 100.0*D;
	
	gl_FragData[0] = vec4(Xp, 1.0);
	gl_FragData[1] = vec4(D, 1.0);
	gl_FragData[2] = state;
	gl_FragData[3] = rgbLambda;
}

