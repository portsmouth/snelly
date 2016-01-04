

precision highp float;

uniform sampler2D Frame;
uniform float Exposure;
varying vec2 vTexCoord;

void main() 
{
	gl_FragColor = texture2D(Frame, vTexCoord); //pow(texture2D(Frame, vTexCoord).rgb*Exposure, vec3(1.0/2.2)), 1.0);
}
