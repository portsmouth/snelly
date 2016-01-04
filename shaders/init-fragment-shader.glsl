
/////////////////////////////////////////////////
// Init fragment shader
/////////////////////////////////////////////////

#extension GL_EXT_draw_buffers : require
precision highp float;

uniform sampler2D RngData;
uniform vec3 EmitterPos;
uniform vec3 EmitterDir;

varying vec2 vTexCoord;

float rand(inout vec4 state)
{
	const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
	const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
	const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
	const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

	vec4 beta = floor(state/q);
	vec4 p = a*(state - beta*q) - beta*r;
	beta = (1.0 - sign(p))*0.5*m;
	state = p + beta;
	return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));
}

void main()
{
	vec4 state = texture2D(RngData, vTexCoord);
	//float rng = rand(state);
	float lambda = 300.0; // + rng*(600.0-300.0);

	vec3 dir = EmitterDir;
	vec3 pos = EmitterPos;
	vec3 rgb = vec3(1.0, 1.0, 1.0);

	gl_FragData[0] = vec4(pos, 1.0);
	gl_FragData[1] = vec4(dir, 1.0);
	gl_FragData[2] = state;
	gl_FragData[3] = vec4(rgb, lambda);
}
