
uniform sampler2D RngData;
uniform sampler2D Spectrum;

uniform vec3 EmitterPos;
uniform vec3 EmitterDir;
uniform float EmitterRadius;
uniform float EmitterSpread; // in degrees

varying vec2 vTexCoord;

void main()
{
	vec4 state = texture2D(RngData, vTexCoord);

	float lambda = 360.0 + (750.0 - 360.0)*vTexCoord.x;
	vec3 rgb = texture2D(Spectrum, vec2(vTexCoord.x, 0.5)).rgb;

	// @todo: make cross-section circular
	vec3 pos = EmitterPos + EmitterRadius*(-vec3(0.5) + vec3(rand(state), rand(state), rand(state)));

	// @todo: ncorrect, do properly (e.g. 90 degree spread != hemisphere with this)
	float spreadAngle = EmitterSpread*(M_PI/360.0);
	vec3 dir = normalize(EmitterDir + spreadAngle*(-vec3(0.5) + vec3(rand(state), rand(state), rand(state))));
	
	gl_FragData[0] = vec4(pos, 1.0);
	gl_FragData[1] = vec4(dir, 1.0);
	gl_FragData[2] = state;
	gl_FragData[3] = vec4(rgb, lambda);
}
