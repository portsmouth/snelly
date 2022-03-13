precision highp float;

uniform sampler2D Radiance;
in vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float contrast;
uniform float saturation;
uniform float hueShift;

out vec4 g_outputColor;

float toneMap(float L)
{
  return L / (1.0 + L);
}

vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));
    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{
    vec3 L = texture(Radiance, vTexCoord).rgb;
    float X = L.x;
    float Y = L.y;
    float Z = L.z;
    
    // convert XYZ tristimulus to linear RGB color space
    vec3 RGB;
    RGB.r =  3.2406*X - 1.5372*Y - 0.4986*Z;
    RGB.g = -0.9689*X + 1.8758*Y + 0.0415*Z;
    RGB.b =  0.0557*X - 0.2040*Y + 1.0570*Z;
    
    // deal with out-of-gamut RGB.
    float delta = -min(0.0, min(min(RGB.r, RGB.g), RGB.b));
    RGB.r += delta;
    RGB.g += delta;
    RGB.b += delta;

    // apply gamma correction to convert linear RGB to sRGB
    RGB = pow(RGB, vec3(invGamma));

    // apply tonemapping
    RGB *= pow(2.0, exposure);
    float R = RGB.r;
    float G = RGB.g;
    float B = RGB.b;
    R = toneMap(R);
    G = toneMap(G);
    B = toneMap(B);



    // apply saturation
    float mean = (R + G + B)/3.0;
    float dR = R - mean;
    float dG = G - mean;
    float dB = B - mean;
    R = mean + sign(dR)*pow(abs(dR), 1.0/saturation);
    G = mean + sign(dG)*pow(abs(dG), 1.0/saturation);
    B = mean + sign(dB)*pow(abs(dB), 1.0/saturation);

    // apply contrast
    dR = R - 0.5;
    dG = G - 0.5;
    dB = B - 0.5;
    R = 0.5 + sign(dR)*pow(abs(dR), 1.0/contrast);
    G = 0.5 + sign(dG)*pow(abs(dG), 1.0/contrast);
    B = 0.5 + sign(dB)*pow(abs(dB), 1.0/contrast);

    vec3 C =vec3(R,G,B);

    // apply hue-shift
    if (hueShift > 0.0)
    {
        vec3 hsv = rgb2hsv(C);
        hsv.r = mod(hsv.r + hueShift, 1.0);
        C = hsv2rgb(hsv);
    }

    g_outputColor = vec4(C, 1.0);
}
