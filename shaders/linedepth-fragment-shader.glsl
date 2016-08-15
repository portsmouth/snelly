
uniform sampler2D PosDataA;
uniform sampler2D PosDataB;
uniform sampler2D Depth;

uniform float camNear;
uniform float camFar;

varying vec2 vTexCoord;
varying float eye_z;


void main() 
{
	float destClipDepth = unpack_depth( texture2D(Depth, vTexCoord) );
	
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






