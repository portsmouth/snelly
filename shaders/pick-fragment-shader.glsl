precision highp float;

layout(location = 0) out vec4 gbuf_pick;

uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees 
uniform float camAspect;
uniform float camAperture;
uniform float camFocalDistance;

uniform float minLengthScale;
uniform float maxLengthScale;
uniform bool maxStepsIsMiss;
uniform vec2 mousePick;

#define MAT_INVAL  -1
#define MAT_DIELE  0
#define MAT_METAL  1
#define MAT_SURFA  2
#define MAT_VOLUM  3

#define M_PI 3.1415926535897932384626433832795

/// GLSL floating point pseudorandom number generator, from
/// "Implementing a Photorealistic Rendering System using GLSL", Toshiya Hachisuka
/// http://arxiv.org/pdf/1505.06022.pdf
float rand(inout vec4 rnd) 
{
    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);
    vec4 beta = floor(rnd/q);
    vec4 p = a*(rnd - beta*q) - beta*r;
    beta = (1.0 - sign(p))*0.5*m;
    rnd = p + beta;
    return fract(dot(rnd/m, vec4(1.0, -1.0, 1.0, -1.0)));
}


//////////////////////////////////////////////////////////////
// Dynamically injected code
//////////////////////////////////////////////////////////////

__DEFINES__

__SHADER__

///////////////////////////////////////////////////////////////////////////////////
// SDF raymarcher
///////////////////////////////////////////////////////////////////////////////////

// find first hit over specified segment
bool traceDistance(in vec3 start, in vec3 dir, float maxDist,
                   inout vec3 hit, inout int material)
{
    float minMarchDist = minLengthScale;

    const float HUGE_VAL = 1.0e20;
    float sdf = HUGE_VAL;
#ifdef HAS_SURFACE
    float sdf_surfa = abs(SDF_SURFACE(start));    sdf = min(sdf, sdf_surfa);
#endif
#ifdef HAS_METAL
    float sdf_metal = abs(SDF_METAL(start));      sdf = min(sdf, sdf_metal);
#endif
#ifdef HAS_DIELECTRIC
    float sdf_diele = abs(SDF_DIELECTRIC(start)); sdf = min(sdf, sdf_diele);
#endif
#ifdef HAS_VOLUME
    float sdf_volum = abs(SDF_VOLUME(start)); sdf = min(sdf, sdf_volum);
#endif

    float InitialSign = sign(sdf);
    float t = 0.0;
    int iters=0;
    for (int n=0; n<__MAX_MARCH_STEPS__; n++)
    {
        // With this formula, the ray advances whether sdf is initially negative or positive --
        // but on crossing the zero isosurface, sdf flips allowing bracketing of the root.
        t += InitialSign * sdf;
        if (t>=maxDist) break;
        vec3 pW = start + t*dir;
        sdf = HUGE_VAL;
#ifdef HAS_SURFACE
        sdf_surfa = abs(SDF_SURFACE(pW));    if (sdf_surfa<minMarchDist) { material = MAT_SURFA; break; } sdf = min(sdf, sdf_surfa);
#endif
#ifdef HAS_METAL
        sdf_metal = abs(SDF_METAL(pW));      if (sdf_metal<minMarchDist) { material = MAT_METAL; break; } sdf = min(sdf, sdf_metal);
#endif
#ifdef HAS_DIELECTRIC
        sdf_diele = abs(SDF_DIELECTRIC(pW)); if (sdf_diele<minMarchDist) { material = MAT_DIELE; break; } sdf = min(sdf, sdf_diele); 
#endif
#ifdef HAS_VOLUME
        sdf_volum = abs(SDF_VOLUME(pW)); if (sdf_volum<minMarchDist) { material = MAT_VOLUM; break; } sdf = min(sdf, sdf_volum); 
#endif
        iters++;
    }
    hit = start + t*dir;
    if (t>=maxDist) return false;
    if (maxStepsIsMiss && iters>=__MAX_MARCH_STEPS__) return false;
    return true;
}

// find first hit along infinite ray
bool traceRay(in vec3 start, in vec3 dir, 
              inout vec3 hit, inout int material, float maxMarchDist)
{
    material = MAT_INVAL;
    return traceDistance(start, dir, maxMarchDist, hit, material);
}

////////////////////////////////////////////////////////////////////////////////
// Picking program: computes first hit distance based on mouse location
////////////////////////////////////////////////////////////////////////////////

void constructPickRay(in vec2 pixel,
                      inout vec3 primaryStart, inout vec3 primaryDir)
{
    // Compute world ray direction for given fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
    float fh = tan(0.5*radians(camFovy)); // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    primaryDir = normalize(camDir + s);
    primaryStart = camPos;
}

void main()
{
    INIT();

    vec2 pixel = mousePick;
    vec3 primaryStart, primaryDir;
    constructPickRay(pixel, primaryStart, primaryDir);

     // Raycast to first hit point
    vec3 pW;
    int hitMaterial;
    bool hit = traceRay(primaryStart, primaryDir, pW, hitMaterial, maxLengthScale);
    float d = maxLengthScale;
    if (hit)
    {
        d = length(pW - primaryStart);
    }

    gbuf_pick = vec4(d, hitMaterial, 0, 0);
}


