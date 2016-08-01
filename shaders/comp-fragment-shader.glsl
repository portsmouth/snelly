
uniform sampler2D Frame;
uniform float Exposure;
varying vec2 vTexCoord;

void main() 
{
	// @todo: expose gamma here in UI.
	// @todo: 'proper' tonemapping here.
	gl_FragColor = vec4(pow(texture2D(Frame, vTexCoord).rgb*Exposure, vec3(1.0/2.2)), 1.0);
}
