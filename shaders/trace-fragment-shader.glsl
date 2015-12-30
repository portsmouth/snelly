
/////////////////////////////////////////////////
// Raytrace fragment shader
/////////////////////////////////////////////////

//#extension GL_EXT_draw_buffers : require
precision mediump float;

//uniform sampler2D PosData;
//uniform sampler2D DirData;
//uniform sampler2D RngData;
//uniform sampler2D RgbData;

void main() 
{
	float i = float(gl_FragCoord.x) / 512.0;
	float j = float(gl_FragCoord.y) / 512.0;

	// Read X, W, F
	// Do:  X' = trace(X, W) given F
	// Write X' into fragment (of ray buffer)
	vec3 Xp = vec3(i, j, 1.0);
	
	// Will write Xp into a texture (via render-to-texture)
	gl_FragColor = vec4(Xp, 1.0);
}

