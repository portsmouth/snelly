
var Pathtracer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;

	// Create a quad VBO for rendering textures
	this.quadVbo = this.createQuadVbo();
	
	this.fbo = new GLU.RenderTarget();
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	var render_canvas = document.getElementById('render-canvas');

	this.shaderSources = GLU.resolveShaderSource(["pathtracer"]);
	this.shaderSources = GLU.resolveShaderSource(["tonemapper"]);

	this.compileShaders();

}

Pathtracer.prototype.createQuadVbo = function()
{
	var vbo = new GLU.VertexBuffer();
	vbo.addAttribute("Position", 3, this.gl.FLOAT, false);
	vbo.addAttribute("TexCoord", 2, this.gl.FLOAT, false);
	vbo.init(4);
	vbo.copy(new Float32Array([
		 1.0,  1.0, 0.0, 1.0, 1.0,
		-1.0,  1.0, 0.0, 0.0, 1.0,
		-1.0, -1.0, 0.0, 0.0, 0.0,
		 1.0, -1.0, 0.0, 1.0, 0.0
	]));
	return vbo;
}

Pathtracer.prototype.reset = function()
{
	this.compileShaders();

	// clear radiance buffer
	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.radianceBuffer, 0);
	this.gl.clear(this.gl.COLOR_BUFFER_BIT);
	this.fbo.unbind();
}


Pathtracer.prototype.compileShaders = function()
{
	console.log('Pathtracer.prototype.compileShaders');

	// Inject code for the current scene and material:
	var sdfCode    = renderer.sceneObj.sdf();

	// Copy the current scene and material routines into the source code
	// of the trace fragment shader
	replacements = {};
	replacements.SDF_FUNC    = sdfCode;
	replacements.MAX_MARCH_STEPS = renderer.maxMarchSteps;

	// shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.pathtraceProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
	this.tonemapProgram =   new GLU.Shader('tonemap', this.shaderSources, null);
}


Pathtracer.prototype.render = function()
{
	console.log('Pathtracer.prototype.render');

	gl.viewport(0, 0, this.width, this.height);

	this.pathtraceProgram.uniform2Fv("resolution", [this.width, this.height]);

	// set camera info 	
	var camPos = renderer.camera.position;
	var camDir = renderer.camera.getWorldDirection();
	var camUp  = renderer.camera.up.transformDirection( this.camera.matrix );
	var camX = camUp.cross(camDir);
	this.pathtraceProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	this.pathtraceProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	this.pathtraceProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	this.pathtraceProgram.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
	this.pathtraceProgram.uniformF("camNear", camera.near);
	this.pathtraceProgram.uniformF("camFar", camera.far);
	this.pathtraceProgram.uniformF("camFovy", camera.fov);

	// trace a path per pixel, progressively updating the HDR radiance into the 'radiance buffer' texture
	this.fbo.drawBuffers(2);
	this.fbo.attachTexture(this.radianceBuffer, 0);
	this.fbo.attachTexture(this.rngTex, 1);
	if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) 
	{
		GLU.fail("[pathtracer] Invalid framebuffer");
	}
	this.radianceTex.bind(0);
	this.rngTex.bind(1);
	this.pathtraceProgram.uniformTexture("Radiance", this.radianceTex);
	this.pathtraceProgram.uniformTexture("RngData", this.rngTex);

	this.pathtraceProgram.bind();
	this.quadVbo.bind();
	this.quadVbo.draw(this.pathtraceProgram, gl.TRIANGLE_FAN);
	this.fbo.unbind();


	// tonemap the radiance buffer to produce the final pixels
	this.fbo.unbind();
	this.radianceTex.bind(0);
	this.tonemapProgram.bind();
	this.tonemapProgram.uniformTexture("Radiance", this.radianceTex);
	this.quadVbo.bind();
	this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);
}



Pathtracer.prototype.resize = function(width, height)
{
	this.width = width;
	this.height = height;

	// Incrementally updated radiance buffer
	this.radianceTex = new GLU.Texture(this.width, this.height, 4, true, false, true, null);

	// Random number seed buffer
	var rngData = new Float32Array(this.width*this.height*4);
	for (var i = 0; i<this.width*this.height; ++i)
	{
		for (var t = 0; t<4; ++t)
		{
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
	}
	this.rngTex = new GLU.Texture(this.width, this.height, 4, true, false, true, rngData);
}


