
/////////////////////////////////////////////////
// Raytrace fragment shader
/////////////////////////////////////////////////

#extension GL_EXT_draw_buffers : require
precision highp float;

uniform sampler2D PosData;
uniform sampler2D DirData;
uniform sampler2D RngData;
uniform sampler2D RgbData;

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


///////////////////////////////
// Distance field
///////////////////////////////

// For now, assume this defines a solid body, whose interior
// is defined by the points with SDF<0.0, with a constant refractive index.


float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) +
         length(max(d,0.0));
}

float sdSphere( vec3 p, vec3 c, float s )
{
  return length(p-c) - s;
}

float sdTorus( vec3 p, float r, float R )
{
  vec2 q = vec2(length(p.xz) - R, p.y);
  return length(q) - r;
}

/*
float opRep( vec3 p, vec3 c )
{
    vec3 q = mod(p,c) - 0.5*c;
    return sdSphere( q, 1.0 );
}
*/

const float inf = 1.0e6;

float sdCross( in vec3 p )
{
  float da = sdBox(p.xyz,vec3(inf,1.0,1.0));
  float db = sdBox(p.yzx,vec3(1.0,inf,1.0));
  float dc = sdBox(p.zxy,vec3(1.0,1.0,inf));
  return min(da,min(db,dc));
}

vec3 map( in vec3 p )
{
   float d = sdBox(p,vec3(1.0));

   float s = 1.0;
   for( int m=0; m<1; m++ )
   {
      vec3 a = mod( p*s, 2.0 )-1.0;
      s *= 2.0;
      vec3 r = abs(1.0 - 3.0*abs(a));

      float da = max(r.x,r.y);
      float db = max(r.y,r.z);
      float dc = max(r.z,r.x);
      float c = (min(da,min(db,dc))-1.0)/s;

      d = max(d,c);
   }

   return vec3(d,1.0,1.0);
}

float DF(vec3 X)
{
	//return map(X/4.0).x;
	float SDF1 = sdTorus(X, 0.8, 1.8);
	float SDF2 = sdSphere(X, vec3(0.0, 0.0, 0.0), 1.99);
	return min(SDF1, SDF2);
}




#define normalEpsilon 5.0e-4

vec3 calcNormal( in vec3 pos )
{
	// Compute normal as gradient of SDF
	vec3 eps = vec3(normalEpsilon, 0.0, 0.0);
	vec3 nor = vec3(
	    DF(pos+eps.xyy) - DF(pos-eps.xyy),
	    DF(pos+eps.yxy) - DF(pos-eps.yxy),
	    DF(pos+eps.yyx) - DF(pos-eps.yyx) );
	return normalize(nor);
}


float sellmeierIor(vec3 b, vec3 c, float lambda) 
{
	// (where lambda is in nanometres)
	float lSq = (lambda*1e-3)*(lambda*1e-3);
	return 1.0 + dot((b*lSq)/(lSq - c), vec3(1.0));
	//return sqrt(2.0);
}

// N is outward normal (from solid to vacuum)
vec3 refract(inout vec3 X, vec3 D, vec3 N, inout vec4 rgbLambda)
{
	float lambda = rgbLambda.w;
	float ior = sqrt(sellmeierIor(vec3(1.0396, 0.2318, 1.0105), 
								  vec3(0.0060, 0.0200, 103.56), 
								  lambda));

	float cosi = dot(-D, N);
	bool entering = cosi > 0.0;
	float ei, et;
	if (entering)
	{
		ei = 1.0; // Incident from vacuum, if entering
		et = ior; // Transmitted to internal medium, if entering
	}
	else
	{
		ei = ior;  // Incident from internal medium, if exiting
		et = 1.0;  // Transmitted to vacuum, if exiting
		N *= -1.0; // Flip normal (so normal is always opposite to incident light direction)
		cosi *= -1.0;
	}

	float sini = sqrt(max(0.0, 1.0 - cosi*cosi));
	float r = ei/et;
	float sint = r * sini;

	// Handle total internal reflection (occurs only if exiting)
	if (sint >= 1.0)
	{
		// Shift X slightly away from surface, towards interior
		X += normalEpsilon*N;
		return D - 2.0*dot(D,N)*N;
	}

	float cost = sqrt(max(0.0, 1.0 - sint*sint));
	rgbLambda.rgb /= r*r;
	
	X -= normalEpsilon*N;
	return r*D + (r*cosi - cost)*N; // transmitted direction
}


#define maxMarchSteps 128
#define minMarchDist 1.0e-4


void raytrace(vec3 X, vec3 D,
			  inout vec3 Xp, inout vec3 Dp,
			  inout vec4 rgbLambda)
{
	float totalDistance = 0.0;
	bool hit = false;
	for (int i=0; i<maxMarchSteps; i++)
	{
		Xp = X + totalDistance*D;
		float dist = abs(DF(Xp));
		totalDistance += dist;
		if (dist < minMarchDist)
		{
			hit = true;
			break;
		}
	}
	if (!hit)
	{
		Dp = D;
	}
	else
	{
		vec3 N = calcNormal(Xp);
		Dp = refract(Xp, D, N, rgbLambda);
	}
}


void main()
{
	vec3 X         = texture2D(PosData, vTexCoord).xyz;
	vec3 D         = texture2D(DirData, vTexCoord).xyz;
	vec4 state     = texture2D(RngData, vTexCoord);
	vec4 rgbLambda = texture2D(RgbData, vTexCoord);

	vec3 Xp, Dp;
	raytrace(X, D, Xp, Dp, rgbLambda);

	gl_FragData[0] = vec4(Xp, 1.0);
	gl_FragData[1] = vec4(Dp, 1.0);
	gl_FragData[2] = state;
	gl_FragData[3] = rgbLambda;
}

