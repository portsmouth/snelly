
uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D Depth;

uniform float camNear;
uniform float camFar;

varying float eye_z;
uniform vec2 resolution;


void main() 
{
	vec2 viewportTexCoords = gl_FragCoord.xy/resolution;
	float destClipDepth = unpack_depth( texture2D(Depth, viewportTexCoords) );
	
	float sourceClipDepth = computeClipDepth(eye_z, camNear, camFar);
	if (sourceClipDepth < destClipDepth)
	{
		gl_FragColor = pack_depth(sourceClipDepth);
	}
	else
	{
		gl_FragColor = pack_depth(destClipDepth);
	}
}






