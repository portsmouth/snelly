

#extension GL_EXT_draw_buffers : require
precision highp float;

uniform sampler2D Radiance;      
uniform vec2 resolution;

#define RADIANCE_TOLERANCE 1.0e-6

float filter(float r2)
{
    return max(0.0, 1.0-r2);
}

void main()
{
    vec2 pixel = gl_FragCoord.xy;

    // Average current radiance over sample set within filter radius
    float maxx = resolution.x-1.0;
    float maxy = resolution.y-1.0;
    vec2 invRes = 1.0/resolution.xy;
    float invWidth = 0.5/float(DOWN_RES);

    vec3 mean = vec3(0.0);
    float norm = 0.0;
    for (int i=-DOWN_RES; i<=DOWN_RES; ++i)
    {
        float _i = max(min(maxx, pixel.x+float(i)), 0.0);
        float u = _i*invRes.x;
        float dx = float(i)*invWidth;
        float dx2 = dx*dx;

        for (int j=-DOWN_RES; j<DOWN_RES; ++j)
        {
            float _j = max(min(maxy, pixel.y+float(j)), 0.0);
            float v = _j*invRes.y;
            float dy = float(j)*invWidth;
            
            float r2 = dx2 + dy*dy;
            float f = filter(r2);
            if (f>1.0e-3)
            {
                vec4 L = texture2D(Radiance, vec2(u,v));  
                float weight = f * L.w;
                mean += weight * L.xyz;
                norm += weight;
            }
        }
    }

    mean /= max(norm, RADIANCE_TOLERANCE);
    gl_FragData[0] = vec4(mean, 1.0);
}








