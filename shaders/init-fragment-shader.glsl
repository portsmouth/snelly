
uniform sampler2D RngData;
uniform sampler2D Spectrum;

uniform vec3 EmitterPos;
uniform vec3 EmitterDir;
uniform float EmitterRadius;
uniform float EmitterSpread; // in degrees

varying vec2 vTexCoord;

void main()
{
	vec4 seed = texture2D(RngData, vTexCoord);

	float lambda = 360.0 + (750.0 - 360.0)*vTexCoord.x;
	vec3 rgb = texture2D(Spectrum, vec2(vTexCoord.x, 0.5)).rgb;

	// Make emission cross-section circular
	float rPos   = EmitterRadius*rand(seed);
	float phiPos = 2.0*M_PI*rand(seed);
	vec3 X = vec3(1.0, 0.0, 0.0);
	vec3 Z = vec3(0.0, 0.0, 1.0);
	vec3 u = cross(Z, EmitterDir);
	if ( length(u) < 1.0e-3 )
	{
		u = cross(X, EmitterDir);
	}
	vec3 v = cross(EmitterDir, u);
	vec3 pos = EmitterPos + rPos*(u*cos(phiPos) + v*sin(phiPos)); 

	// Emit in a cone with the given spread
	float spreadAngle = abs(EmitterSpread)*M_PI/180.0;
	float rDir = min(tan(spreadAngle), 1.0e6);
	float phiDir = 2.0*M_PI*rand(seed);
	vec3 dir = normalize(EmitterDir + rDir*(u*cos(phiDir) + v*sin(phiDir)));
	
	gl_FragData[0] = vec4(pos, 1.0);
	gl_FragData[1] = vec4(dir, 1.0);
	gl_FragData[2] = seed;
	gl_FragData[3] = vec4(rgb, lambda);
}
