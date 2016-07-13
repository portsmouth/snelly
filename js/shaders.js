var Shaders = {

'comp-fragment-shader': 
	'precision highp float;\n' +
	'\n' +
	'uniform sampler2D Frame;\n' +
	'uniform float Exposure;\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_FragColor = vec4(pow(texture2D(Frame, vTexCoord).rgb*Exposure, vec3(1.0/2.2)), 1.0);\n' +
	'}\n',

'comp-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
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
	'/////////////////////////////////////////////////\n' +
	'// Init fragment shader\n' +
	'/////////////////////////////////////////////////\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'precision highp float;\n' +
	'\n' +
	'uniform sampler2D RngData;\n' +
	'uniform sampler2D Spectrum;\n' +
	'\n' +
	'uniform vec3 EmitterPos;\n' +
	'uniform vec3 EmitterDir;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
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
	'\n' +
	'void main()\n' +
	'{\n' +
	'	vec4 state = texture2D(RngData, vTexCoord);\n' +
	'\n' +
	'	float lambda = 360.0 + (750.0 - 360.0)*vTexCoord.x;\n' +
	'	vec3 rgb = texture2D(Spectrum, vec2(vTexCoord.x, 0.5)).rgb;\n' +
	'\n' +
	'	vec3 pos = EmitterPos + 0.5*(-vec3(0.5) + vec3(rand(state), rand(state), rand(state)));\n' +
	'	vec3 dir = normalize(EmitterDir + 0.01*vec3(rand(state), rand(state), rand(state)));\n' +
	'	\n' +
	'	gl_FragData[0] = vec4(pos, 1.0);\n' +
	'	gl_FragData[1] = vec4(dir, 1.0);\n' +
	'	gl_FragData[2] = state;\n' +
	'	gl_FragData[3] = vec4(rgb, lambda);\n' +
	'}\n',

'init-vertex-shader': 
	'/////////////////////////////////////////////////\n' +
	'// Init vertex shader\n' +
	'/////////////////////////////////////////////////\n' +
	'\n' +
	'precision highp float;\n' +
	'\n' +
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
	'/////////////////////////////////////////////////\n' +
	'// Viewport fragment shader\n' +
	'/////////////////////////////////////////////////\n' +
	'\n' +
	'precision highp float;\n' +
	'\n' +
	'varying vec3 vColor;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_FragColor = vec4(vColor, 1.0);\n' +
	'}\n',

'line-vertex-shader': 
	'/////////////////////////////////////////////////\n' +
	'// Line vertex shader\n' +
	'/////////////////////////////////////////////////\n' +
	'\n' +
	'precision highp float;\n' +
	'\n' +
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
	'uniform sampler2D Frame;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
	'\n' +
	'void main() \n' +
	'{\n' +
	'	gl_FragColor = vec4(texture2D(Frame, vTexCoord).rgb, 1.0);\n' +
	'}\n',

'pass-vertex-shader': 
	'precision highp float;\n' +
	'\n' +
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
	'/////////////////////////////////////////////////\n' +
	'// Raytrace fragment shader\n' +
	'/////////////////////////////////////////////////\n' +
	'\n' +
	'#extension GL_EXT_draw_buffers : require\n' +
	'precision highp float;\n' +
	'\n' +
	'uniform sampler2D PosData;\n' +
	'uniform sampler2D DirData;\n' +
	'uniform sampler2D RngData;\n' +
	'uniform sampler2D RgbData;\n' +
	'\n' +
	'varying vec2 vTexCoord;\n' +
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
	'\n' +
	'\n' +
	'///////////////////////////////\n' +
	'// Distance field\n' +
	'///////////////////////////////\n' +
	'\n' +
	'// For now, assume this defines a solid body, whose interior\n' +
	'// is defined by the points with SDF<0.0, with a constant refractive index.\n' +
	'\n' +
	'\n' +
	'float sdBox( vec3 p, vec3 b )\n' +
	'{\n' +
	'  vec3 d = abs(p) - b;\n' +
	'  return min(max(d.x,max(d.y,d.z)),0.0) +\n' +
	'         length(max(d,0.0));\n' +
	'}\n' +
	'\n' +
	'float sdSphere( vec3 p, vec3 c, float s )\n' +
	'{\n' +
	'  return length(p-c) - s;\n' +
	'}\n' +
	'\n' +
	'float sdTorus( vec3 p, float r, float R )\n' +
	'{\n' +
	'  vec2 q = vec2(length(p.xz) - R, p.y);\n' +
	'  return length(q) - r;\n' +
	'}\n' +
	'\n' +
	'/*\n' +
	'float opRep( vec3 p, vec3 c )\n' +
	'{\n' +
	'    vec3 q = mod(p,c) - 0.5*c;\n' +
	'    return sdSphere( q, 1.0 );\n' +
	'}\n' +
	'*/\n' +
	'\n' +
	'const float inf = 1.0e6;\n' +
	'\n' +
	'float sdCross( in vec3 p )\n' +
	'{\n' +
	'  float da = sdBox(p.xyz,vec3(inf,1.0,1.0));\n' +
	'  float db = sdBox(p.yzx,vec3(1.0,inf,1.0));\n' +
	'  float dc = sdBox(p.zxy,vec3(1.0,1.0,inf));\n' +
	'  return min(da,min(db,dc));\n' +
	'}\n' +
	'\n' +
	'vec3 map( in vec3 p )\n' +
	'{\n' +
	'   float d = sdBox(p,vec3(1.0));\n' +
	'\n' +
	'   float s = 1.0;\n' +
	'   for( int m=0; m<1; m++ )\n' +
	'   {\n' +
	'      vec3 a = mod( p*s, 2.0 )-1.0;\n' +
	'      s *= 2.0;\n' +
	'      vec3 r = abs(1.0 - 3.0*abs(a));\n' +
	'\n' +
	'      float da = max(r.x,r.y);\n' +
	'      float db = max(r.y,r.z);\n' +
	'      float dc = max(r.z,r.x);\n' +
	'      float c = (min(da,min(db,dc))-1.0)/s;\n' +
	'\n' +
	'      d = max(d,c);\n' +
	'   }\n' +
	'\n' +
	'   return vec3(d,1.0,1.0);\n' +
	'}\n' +
	'\n' +
	'float DF(vec3 X)\n' +
	'{\n' +
	'	return map(X/4.0).x;\n' +
	'	//float SDF1 = sdTorus(X, 0.8, 1.8);\n' +
	'	//float SDF2 = sdSphere(X, vec3(0.0, 0.0, 0.0), 1.99);\n' +
	'	//return min(SDF1, SDF2);\n' +
	'}\n' +
	'\n' +
	'\n' +
	'\n' +
	'\n' +
	'#define normalEpsilon 5.0e-4\n' +
	'\n' +
	'vec3 calcNormal( in vec3 pos )\n' +
	'{\n' +
	'	// Compute normal as gradient of SDF\n' +
	'	vec3 eps = vec3(normalEpsilon, 0.0, 0.0);\n' +
	'	vec3 nor = vec3(\n' +
	'	    DF(pos+eps.xyy) - DF(pos-eps.xyy),\n' +
	'	    DF(pos+eps.yxy) - DF(pos-eps.yxy),\n' +
	'	    DF(pos+eps.yyx) - DF(pos-eps.yyx) );\n' +
	'	return normalize(nor);\n' +
	'}\n' +
	'\n' +
	'\n' +
	'float sellmeierIor(vec3 b, vec3 c, float lambda) \n' +
	'{\n' +
	'	// (where lambda is in nanometres)\n' +
	'	float lSq = (lambda*1e-3)*(lambda*1e-3);\n' +
	'	return 1.0 + dot((b*lSq)/(lSq - c), vec3(1.0));\n' +
	'	//return sqrt(2.0);\n' +
	'}\n' +
	'\n' +
	'// N is outward normal (from solid to vacuum)\n' +
	'vec3 refract(inout vec3 X, vec3 D, vec3 N, inout vec4 rgbLambda)\n' +
	'{\n' +
	'	float lambda = rgbLambda.w;\n' +
	'	float ior = sqrt(sellmeierIor(vec3(1.0396, 0.2318, 1.0105), \n' +
	'								  vec3(0.0060, 0.0200, 103.56), \n' +
	'								  lambda));\n' +
	'\n' +
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
	'	float sini = sqrt(max(0.0, 1.0 - cosi*cosi));\n' +
	'	float r = ei/et;\n' +
	'	float sint = r * sini;\n' +
	'\n' +
	'	// Handle total internal reflection (occurs only if exiting)\n' +
	'	if (sint >= 1.0)\n' +
	'	{\n' +
	'		// Shift X slightly away from surface, towards interior\n' +
	'		X += normalEpsilon*N;\n' +
	'		return D - 2.0*dot(D,N)*N;\n' +
	'	}\n' +
	'\n' +
	'	float cost = sqrt(max(0.0, 1.0 - sint*sint));\n' +
	'	rgbLambda.rgb /= r*r;\n' +
	'	\n' +
	'	X -= normalEpsilon*N;\n' +
	'	return r*D + (r*cosi - cost)*N; // transmitted direction\n' +
	'}\n' +
	'\n' +
	'\n' +
	'#define maxMarchSteps 128\n' +
	'#define minMarchDist 1.0e-4\n' +
	'\n' +
	'\n' +
	'void raytrace(vec3 X, vec3 D,\n' +
	'			  inout vec3 Xp, inout vec3 Dp,\n' +
	'			  inout vec4 rgbLambda)\n' +
	'{\n' +
	'	float totalDistance = 0.0;\n' +
	'	bool hit = false;\n' +
	'	for (int i=0; i<maxMarchSteps; i++)\n' +
	'	{\n' +
	'		Xp = X + totalDistance*D;\n' +
	'		float dist = abs(DF(Xp));\n' +
	'		totalDistance += dist;\n' +
	'		if (dist < minMarchDist)\n' +
	'		{\n' +
	'			hit = true;\n' +
	'			break;\n' +
	'		}\n' +
	'	}\n' +
	'	if (!hit)\n' +
	'	{\n' +
	'		Dp = D;\n' +
	'	}\n' +
	'	else\n' +
	'	{\n' +
	'		vec3 N = calcNormal(Xp);\n' +
	'		Dp = refract(Xp, D, N, rgbLambda);\n' +
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
	'	vec3 Xp, Dp;\n' +
	'	raytrace(X, D, Xp, Dp, rgbLambda);\n' +
	'\n' +
	'	gl_FragData[0] = vec4(Xp, 1.0);\n' +
	'	gl_FragData[1] = vec4(Dp, 1.0);\n' +
	'	gl_FragData[2] = state;\n' +
	'	gl_FragData[3] = rgbLambda;\n' +
	'}\n',

'trace-vertex-shader': 
	'/////////////////////////////////////////////////\n' +
	'// Raytrace vertex shader\n' +
	'/////////////////////////////////////////////////\n' +
	'\n' +
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