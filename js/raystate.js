

var RayState = function(size) 
{
	var posData = new Float32Array(size*size*4); // ray position
	var dirData = new Float32Array(size*size*4); // ray direction
	var rngData = new Float32Array(size*size*4); // Random number seed
	var rgbData = new Float32Array(size*size*4); // Ray color, and wavelength

	for (var i = 0; i<size*size; ++i)
	{
		dirData[i*4 + 0] = 1.0;
		dirData[i*4 + 1] = 0.0;
		dirData[i*4 + 2] = 0.0;
		dirData[i*4 + 3] = 0.0;
		for (var t = 0; t<4; ++t)
		{
			rgbData[i*4 + t] = Math.random();
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
	}

	this.posTex = new GLU.Texture(size, size, 4, true, false, true, posData);
	this.dirTex = new GLU.Texture(size, size, 4, true, false, true, dirData);
	this.rngTex = new GLU.Texture(size, size, 4, true, false, true, rngData);
	this.rgbTex = new GLU.Texture(size, size, 4, true, false, true, rgbData);
}

RayState.prototype.bind = function(shader)
{
	this.posTex.bind(0);
	this.dirTex.bind(1);
	this.rngTex.bind(2);
	this.rgbTex.bind(3);
	shader.uniformTexture("PosData", this.posTex);
	shader.uniformTexture("DirData", this.dirTex);
	shader.uniformTexture("RngData", this.rngTex);
	shader.uniformTexture("RgbData", this.rgbTex);
}

RayState.prototype.attach = function(fbo)
{
	var gl = GLU.gl;
	fbo.attachTexture(this.posTex, 0);
	fbo.attachTexture(this.dirTex, 1);
	fbo.attachTexture(this.rngTex, 2);
	fbo.attachTexture(this.rgbTex, 3);
	if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) 
	{
		GLU.fail("Invalid framebuffer");
	}
}

RayState.prototype.detach = function(fbo)
{
	var gl = GLU.gl;
	fbo.detachTexture(0);
	fbo.detachTexture(1);
	fbo.detachTexture(2);
	fbo.detachTexture(3);
}


