
uniform sampler2D WaveBuffer;
varying vec2 vTexCoord;

void main() 
{
	gl_FragColor = vec4(texture2D(WaveBuffer, vTexCoord).rgba);
}
