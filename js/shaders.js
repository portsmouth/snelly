var Shaders = {

'comp-fragment-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'uniform sampler2D Frame;\n' +
	'uniform float Exposure;\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	// @todo: expose gamma here in UI\n' +
	'	gl_FragColor = vec4(pow(texture2D(Frame, vTexCoord).rgb*Exposure, vec3(1.0/2.2)), 1.0);\n' +
	'}\n',

'comp-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'attribute vec3 Position;\n' +
	'attribute vec2 TexCoord;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main(void)\n' +
	'{\n' +
	'	gl_Position = vec4(Position, 1.0);\n' +
	'	vTexCoord = TexCoord;\n' +
	'}\n',

'init-fragment-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'uniform sampler2D RngData;\n' +
	'uniform sampler2D Spectrum;\n' +
	'\n' +
	'uniform vec3 EmitterPos;\n' +
	'uniform vec3 EmitterDir;\n' +
	'uniform float EmitterRadius;\n' +
	'uniform float EmitterSpread; // in degrees\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main()\n' +
	'{\n' +
	'	vec4 state = texture2D(RngData, vTexCoord);\n' +
	'\n' +
	'	float lambda = 360.0 + (750.0 - 360.0)*vTexCoord.x;\n' +
	'	vec3 rgb = texture2D(Spectrum, vec2(vTexCoord.x, 0.5)).rgb;\n' +
	'\n' +
	'	// @todo: make cross-section circular\n' +
	'	vec3 pos = EmitterPos + EmitterRadius*(-vec3(0.5) + vec3(rand(state), rand(state), rand(state)));\n' +
	'\n' +
	'	// @todo: ncorrect, do properly (e.g. 90 degree spread != hemisphere with this)\n' +
	'	float spreadAngle = EmitterSpread*(M_PI/360.0);\n' +
	'	vec3 dir = normalize(EmitterDir + spreadAngle*(-vec3(0.5) + vec3(rand(state), rand(state), rand(state))));\n' +
	'	\n' +
	'	gl_FragData[0] = vec4(pos, 1.0);\n' +
	'	gl_FragData[1] = vec4(dir, 1.0);\n' +
	'	gl_FragData[2] = state;\n' +
	'	gl_FragData[3] = vec4(rgb, lambda);\n' +
	'}\n',

'init-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'attribute vec3 Position;\n' +
	'attribute vec2 TexCoord;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_Position = vec4(Position, 1.0);\n' +
	'	vTexCoord = TexCoord;\n' +
	'}\n',

'line-fragment-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'varying vec3 vColor;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_FragColor = vec4(vColor, 1.0);\n' +
	'}\n',

'line-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'uniform sampler2D PosDataA;\n' +
	'uniform sampler2D PosDataB;\n' +
	'uniform sampler2D RgbData;\n' +
	'\n' +
	'uniform mat4 u_projectionMatrix;\n' +
	'uniform mat4 u_modelViewMatrix;\n' +
	'\n' +
	'attribute vec3 TexCoord;\n' +
	'varying vec3 vColor;\n' +
	'\n' +
	'void main()\n' +
	'{\n' +
	'	// Textures A and B contain line segment start and end points respectively\n' +
	'	vec3 posA = texture2D(PosDataA, TexCoord.xy).xyz;\n' +
	'	vec3 posB = texture2D(PosDataB, TexCoord.xy).xyz;\n' +
	'\n' +
	'	// Line segment vertex position\n' +
	'	vec3 pos = mix(posA, posB, TexCoord.z);\n' +
	'\n' +
	'	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(pos, 1.0);\n' +
	'	vColor = texture2D(RgbData, TexCoord.xy).rgb;\n' +
	'}\n',

'pass-fragment-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'uniform sampler2D Frame;\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_FragColor = vec4(texture2D(Frame, vTexCoord).rgb, 1.0);\n' +
	'}\n',

'pass-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'attribute vec3 Position;\n' +
	'attribute vec2 TexCoord;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main(void)\n' +
	'{\n' +
	'	gl_Position = vec4(Position, 1.0);\n' +
	'	vTexCoord = TexCoord;\n' +
	'}\n',

'trace-fragment-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'uniform sampler2D PosData;\n' +
	'uniform sampler2D DirData;\n' +
	'uniform sampler2D RngData;\n' +
	'uniform sampler2D RgbData;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'// @todo: epsilons should be relative to a scene scale\n' +
	'#define normalEpsilon 5.0e-4\n' +
	'\n' +
	'// @todo: should be user params\n' +
	'#define maxMarchSteps 256\n' +
	'#define minMarchDist 1.0e-4\n' +
	'#define maxDist 1.0e6\n' +
	'\n' +
	'\n' +
	'// N is outward normal (from medium to vacuum)\n' +
	'// Returns radiance gain on reflection (i.e., 1)\n' +
	'float Reflect(inout vec3 X, inout vec3 D, vec3 N)\n' +
	'{\n' +
	'	float cosi = dot(-D, N);\n' +
	'	bool entering = cosi > 0.0;\n' +
	'	if (!entering)\n' +
	'	{\n' +
	'		N *= -1.0; // Flip normal (so normal is always opposite to incident light direction)\n' +
	'	}\n' +
	'	// Reflect direction about normal, and displace ray start into reflected halfspace:\n' +
	'	X += normalEpsilon*N;\n' +
	'	D -= 2.0*N*dot(D,N);\n' +
	'	return 1.0;\n' +
	'}\n' +
	'\n' +
	'\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'// Dielectric formulae\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'\n' +
	'\n' +
	'// N is outward normal (from medium to vacuum).\n' +
	'// Returns radiance gain on transmission\n' +
	'float Transmit(inout vec3 X, inout vec3 D, vec3 N, float ior)\n' +
	'{\n' +
	'	// This applies of course only to dielectrics\n' +
	'	float cosi = dot(-D, N);\n' +
	'	bool entering = cosi > 0.0;\n' +
	'	float ei, et;\n' +
	'	if (entering)\n' +
	'	{\n' +
	'		ei = 1.0; // Incident from vacuum, if entering\n' +
	'		et = ior; // Transmitted to internal medium, if entering\n' +
	'	}\n' +
	'	else\n' +
	'	{\n' +
	'		ei = ior;  // Incident from internal medium, if exiting\n' +
	'		et = 1.0;  // Transmitted to vacuum, if exiting\n' +
	'		N *= -1.0; // Flip normal (so normal is always opposite to incident light direction)\n' +
	'		cosi *= -1.0;\n' +
	'	}\n' +
	'\n' +
	'	float r = ei/et;\n' +
	'\n' +
	'	// Compute sint from Snells law:\n' +
	'	float sint = r * sqrt(max(0.0, 1.0 - cosi*cosi));\n' +
	'\n' +
	'	// sint<=1.0 guaranteed as total internal reflection already handled\n' +
	'	float cost = sqrt(max(0.0, 1.0 - sint*sint)); \n' +
	'\n' +
	'	// Displace ray start into transmitted halfspace\n' +
	'	X -= normalEpsilon*N;\n' +
	'\n' +
	'	// Set transmitted direction\n' +
	'	D = r*D + (r*cosi - cost)*N; \n' +
	'\n' +
	'	// Transmitted radiance gets scaled by the square of the ratio of transmitted to incident IOR:\n' +
	'	return 1.0 / (r*r);\n' +
	'}\n' +
	'\n' +
	'\n' +
	'// D is direction of incident light\n' +
	'// N is normal pointing from medium to vacuum\n' +
	'float reflectionDielectric(vec3 D, vec3 N, float ior)\n' +
	'{\n' +
	'	float cosi = dot(-D, N);\n' +
	'	bool entering = cosi > 0.0;\n' +
	'	float ei, et;\n' +
	'	if (entering)\n' +
	'	{\n' +
	'		ei = 1.0; // Incident from vacuum, if entering\n' +
	'		et = ior; // Transmitted to internal medium, if entering\n' +
	'	}\n' +
	'	else\n' +
	'	{\n' +
	'		ei = ior;  // Incident from internal medium, if exiting\n' +
	'		et = 1.0;  // Transmitted to vacuum, if exiting\n' +
	'		N *= -1.0; // Flip normal (so normal is always opposite to incident light direction)\n' +
	'		cosi *= -1.0;\n' +
	'	}\n' +
	'\n' +
	'	// Compute sint from Snells law:\n' +
	'	float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));\n' +
	'\n' +
	'	// Handle total internal reflection\n' +
	'	if (sint >= 1.0) return 1.0;\n' +
	'	\n' +
	'	float cost = sqrt(max(0.0, 1.0 - sint*sint));\n' +
	'	float cosip = abs(cosi);\n' +
	'	float rParallel      = ( et*cosip - ei*cost) / ( et*cosip + ei*cost);\n' +
	'	float rPerpendicular = ( ei*cosip - et*cost) / ( ei*cosip + et*cost);\n' +
	'	return 0.5 * (rParallel*rParallel + rPerpendicular*rPerpendicular);\n' +
	'}\n' +
	'\n' +
	'\n' +
	'// D is direction of incident light\n' +
	'// N is normal pointing from medium to vacuum\n' +
	'// returns radiance multiplier (1 for reflection, varying for transmission)\n' +
	'float sampleDielectric(inout vec3 X, inout vec3 D, vec3 N, float ior)\n' +
	'{\n' +
	'	float R = reflectionDielectric(D, N, ior);\n' +
	'	float rnd = rand(state);\n' +
	'	if (R >= rnd) // make reflectProb = R\n' +
	'	{\n' +
	'		// we must multiply subsequent radiance by factor (1.0 / reflectProb) to get correct counting\n' +
	'		// but as we chose reflectProb = R, this cancels with R, so that on reflection we should leave the radiance unchanged. \n' +
	'		return Reflect(X, D, N);\n' +
	'	}\n' +
	'	else // refraction, prob = (1-R)\n' +
	'	{\n' +
	'		// we must multiply subsequent radiance by factor (1.0 / transmitProb) to get correct counting\n' +
	'		// but as we chose transmitProb = 1 - reflectProb = 1 - R, this cancels with (1-R) in the numerator, so that\n' +
	'		// on transmission, we should just multiply the radiance by the (et/ei)^2 gain factor (done inside transmission function)\n' +
	'		return Transmit(X, D, N, ior, gain);\n' +
	'	}\n' +
	'}\n' +
	'\n' +
	'\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'// Metal formulae\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'\n' +
	'// D is direction of incident light\n' +
	'// N is normal pointing from medium to vacuum\n' +
	'float reflectionMetal(vec3 D, vec3 N, float ior, float k)\n' +
	'{\n' +
	'	float cosi = dot(-D, N);\n' +
	'	float cosip = abs(cosi);\n' +
	'	float cosi2 = cosip * cosip;\n' +
	'	float tmp = (ior*ior + k*k) * cosi2;\n' +
	'	float twoEtaCosi = 2.0 * ior * cosip;\n' +
	'	float Rparl2 = (tmp - twoEtaCosi + 1.0) / (tmp + twoEtaCosi + 1.0);\n' +
	'	float tmp_f = ior*ior + k*k;\n' +
	'	float Rperp2 = (tmp_f - twoEtaCosi + cosi2) / (tmp_f + twoEtaCosi + cosi2);\n' +
	'	return 0.5 * (Rparl2 + Rperp2);\n' +
	'}\n' +
	'\n' +
	'// D is direction of incident light\n' +
	'// N is normal pointing from medium to vacuum\n' +
	'// returns radiance multiplier (varying for reflection, 0 for 'transmission' i.e. absorption)\n' +
	'float sampleMetal(inout vec3 X, inout vec3 D, vec3 N, float ior, float k)\n' +
	'{\n' +
	'	float R = reflectionMetal(D, N, ior, k);\n' +
	'	if (R >= rnd) // make reflectProb = R\n' +
	'	{\n' +
	'		// we must multiply subsequent radiance by factor (1.0 / reflectProb) to get correct counting\n' +
	'		// but as we chose reflectProb = R, this cancels with R, so that on reflection we should leave the radiance unchanged. \n' +
	'		return Reflect(X, D, N);\n' +
	'	}\n' +
	'	else // absorption prob = (1-R)\n' +
	'	{\n' +
	'		return 0.0;\n' +
	'	}\n' +
	'}\n' +
	'\n' +
	'\n' +
	'\n' +
	'\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'// @todo: paste this code in dynamically, based on current scene\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'\n' +
	'uniform float _radius;                \n' +
	'float SDF(vec3 X)                     \n' +
	'{                                     \n' +
	'	return length(X) - _radius;       \n' +
	'}                                     \n' +
	'\n' +
	'\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'// @todo: paste this code in dynamically, based on current material\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'\n' +
	'float ior(float lnm) \n' +
	'{\n' +
	'	return 1.5;\n' +
	'}\n' +
	'\n' +
	'float sample(inout vec3 X, inout vec3 D, vec3 N, float lnm)\n' +
	'{\n' +
	'	return sampleDielectric(X, D, N, ior(lnm));\n' +
	'}\n' +
	'\n' +
	'\n' +
	'\n' +
	'\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'// Main SDF tracing loop\n' +
	'//////////////////////////////////////////////////////////////\n' +
	'\n' +
	'\n' +
	'vec3 calcNormal(in vec3 X)\n' +
	'{\n' +
	'	// Compute normal as gradient of SDF\n' +
	'	float eps = normalEpsilon;\n' +
	'	vec3 N = vec3( SDF(X+eps) - SDF(X-eps),\n' +
	'				   SDF(X+eps) - SDF(X-eps),\n' +
	'				   SDF(X+eps) - SDF(X-eps) );\n' +
	'	return normalize(N);\n' +
	'}\n' +
	'\n' +
	'\n' +
	'void raytrace(inout vec3 X, inout vec3 D,\n' +
	'			  inout vec4 rgbLambda, \n' +
	'			  inout vec4 state)\n' +
	'{\n' +
	'	const float radianceEpsilon = 1.0e-7;\n' +
	'	if ( length(rgbLambda.rgb) < radianceEpsilon ) return;\n' +
	'\n' +
	'	float totalDist = 0.0;\n' +
	'	bool hit = false;\n' +
	'	for (int i=0; i<maxMarchSteps; i++)\n' +
	'	{\n' +
	'		X += totalDist*D;\n' +
	'		float dist = abs( SDF(X) );\n' +
	'		totalDist += dist;\n' +
	'		if (dist < minMarchDist)\n' +
	'		{\n' +
	'			hit = true;\n' +
	'			break;\n' +
	'		}\n' +
	'	}\n' +
	'\n' +
	'	if (!hit)\n' +
	'	{\n' +
	'		X += maxDist*D;\n' +
	'		rgbLambda.rgb *= 0.0; // terminate ray\n' +
	'	}\n' +
	'	else\n' +
	'	{\n' +
	'		// Hit the surface. Calculate normal there:\n' +
	'		vec3 N = calcNormal(Xp);\n' +
	'\n' +
	'		// Sample new direction, and update radiance accordingly:\n' +
	'		float lambda = rgbLambda.w;\n' +
	'		rgbLambda.rgb *= sample(X, D, N, lambda);\n' +
	'	}\n' +
	'}\n' +
	'\n' +
	'\n' +
	'void main()\n' +
	'{\n' +
	'	vec3 X         = texture2D(PosData, vTexCoord).xyz;\n' +
	'	vec3 D         = texture2D(DirData, vTexCoord).xyz;\n' +
	'	vec4 state     = texture2D(RngData, vTexCoord);\n' +
	'	vec4 rgbLambda = texture2D(RgbData, vTexCoord);\n' +
	'\n' +
	'	raytrace(X, D, rgbLambda, state);\n' +
	'\n' +
	'	gl_FragData[0] = vec4(X, 1.0);\n' +
	'	gl_FragData[1] = vec4(D, 1.0);\n' +
	'	gl_FragData[2] = state;\n' +
	'	gl_FragData[3] = rgbLambda;\n' +
	'}\n',

'trace-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'\n' +
	'#define M_PI 3.1415926535897932384626433832795\n' +
	'\n' +
	'\n' +
	'float rand(inout vec4 state) \n' +
	'{\n' +
	'    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);\n' +
	'    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);\n' +
	'    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);\n' +
	'    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);\n' +
	'\n' +
	'    vec4 beta = floor(state/q);\n' +
	'    vec4 p = a*(state - beta*q) - beta*r;\n' +
	'    beta = (1.0 - sign(p))*0.5*m;\n' +
	'    state = p + beta;\n' +
	'    return fract(dot(state/m, vec4(1.0, -1.0, 1.0, -1.0)));\n' +
	'}\n' +
	'attribute vec3 Position;\n' +
	'attribute vec2 TexCoord;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_Position = vec4(Position, 1.0);\n' +
	'	vTexCoord = TexCoord;\n' +
	'}\n',

}