

var SurfaceRendererState = function(width, height)
{
	var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count 
	var rngData      = new Float32Array(width*height*4); // Random number seed
	for (var i = 0; i<width*height; ++i)
	{
		for (var t = 0; t<4; ++t)
		{
			rngData[i*4 + t] = Math.random();//*4194167.0;
		}
	}
	this.radianceTex = new GLU.Texture(width, height, 4, true, false, true, radianceData);
	this.rngTex      = new GLU.Texture(width, height, 4, true, false, true, rngData);
}

SurfaceRendererState.prototype.bind = function(shader)
{
	this.radianceTex.bind(0);
	this.rngTex.bind(1);
	shader.uniformTexture("Radiance", this.radianceTex);
	shader.uniformTexture("RngData", this.rngTex);
}

SurfaceRendererState.prototype.attach = function(fbo)
{
	var gl = GLU.gl;
	fbo.attachTexture(this.radianceTex, 0);
	fbo.attachTexture(this.rngTex, 1);
	if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) 
	{
		GLU.fail("Invalid framebuffer");
	}
}

SurfaceRendererState.prototype.detach = function(fbo)
{
	var gl = GLU.gl;
	fbo.detachTexture(0);
	fbo.detachTexture(1);
}

SurfaceRendererState.prototype.clear = function(fbo)
{
	// clear radiance buffer
	var gl = GLU.gl;
	fbo.bind();
	fbo.drawBuffers(1);
	fbo.attachTexture(this.radianceTex, 0);
	gl.clear(gl.COLOR_BUFFER_BIT);
	fbo.unbind();
}


var SurfaceRenderer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;

	this.quadVbo = this.createQuadVbo();
	this.fbo = new GLU.RenderTarget();
	
	// Initialize pathtracing textures and fbo:
	this.currentState = 0;
	this.pathStates = [new SurfaceRendererState(this.width, this.height), 
					   new SurfaceRendererState(this.width, this.height)];

	this.maxMarchSteps = 512;

	// Load shaders
	this.shaderSources = GLU.resolveShaderSource(["pathtracer", "tonemapper"]);
	this.compileShaders();

	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	// Trigger initial buffer generation
	this.resize(this.width, this.height);
}

SurfaceRenderer.prototype.createQuadVbo = function()
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

SurfaceRenderer.prototype.reset = function()
{
	this.compileShaders();
	this.currentState = 0;
	this.pathStates[this.currentState].clear(this.fbo);
}


SurfaceRenderer.prototype.compileShaders = function()
{
	// Inject code for the current scene SDF:
	var sceneObj = snelly.getLoadedScene();
	if (sceneObj == null) return;
	var sdfCode = sceneObj.sdf();

	// Copy the current scene and material routines into the source code
	// of the trace fragment shader
	replacements = {};
	replacements.SDF_FUNC        = sdfCode;
	replacements.MAX_MARCH_STEPS = this.maxMarchSteps;

	// shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.pathtraceProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
	this.tonemapProgram =   new GLU.Shader('tonemapper', this.shaderSources, null);
}


SurfaceRenderer.prototype.render = function()
{
	var sceneObj = snelly.getLoadedScene();
	if (sceneObj == null) return;

	var gl = this.gl;
	gl.viewport(0, 0, this.width, this.height);

	this.pathtraceProgram.bind();

	// set camera info 	
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camUp = camera.up.clone();
	camUp.transformDirection( camera.matrixWorld );

	// Camera update
	camera.near = 1.0e-3*sceneObj.getScale();
	camera.far  = 1.0e3*sceneObj.getScale();
	
	var camX = new THREE.Vector3();
	camX.crossVectors(camUp, camDir);

	this.pathtraceProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	this.pathtraceProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	this.pathtraceProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	this.pathtraceProgram.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
	this.pathtraceProgram.uniformF("camNear", camera.near);
	this.pathtraceProgram.uniformF("camFar", camera.far);
	this.pathtraceProgram.uniformF("camFovy", camera.fov);
	this.pathtraceProgram.uniformF("camZoom", camera.zoom);
	this.pathtraceProgram.uniformF("camAspect", camera.aspect);
	this.pathtraceProgram.uniform2Fv("resolution", [this.width, this.height]);
	this.pathtraceProgram.uniformF("SceneScale", sceneObj.getScale()); 

	var current = this.currentState;
	var next    = 1 - current;

	gl.disable(gl.BLEND);

	this.fbo.bind();
	this.fbo.drawBuffers(2);

	// Read data from the 'current' state
	this.pathStates[current].bind(this.pathtraceProgram);

	// Write data into the 'next' state
	this.pathStates[next].attach(this.fbo);

	sceneObj.syncShader(this.pathtraceProgram); 

	// Trace one path per pixel
	this.quadVbo.bind();
	this.quadVbo.draw(this.pathtraceProgram, gl.TRIANGLE_FAN);
	this.fbo.unbind();

	// Tonemap the 'current' radiance buffer to produce this frame's pixels
	this.tonemapProgram.bind();
	this.tonemapProgram.uniformF("exposure", 1.0);
	this.tonemapProgram.uniformF("invGamma", 1.0/2.2);

	var radianceTexCurrent = this.pathStates[current].radianceTex;
	radianceTexCurrent.bind(0);
	this.tonemapProgram.uniformTexture("Radiance", radianceTexCurrent);
	
	gl.blendFunc(gl.ONE, gl.ONE);
	gl.enable(gl.BLEND);

	this.quadVbo.bind();
	this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);

	this.fbo.unbind();
	gl.disable(gl.BLEND);

	// Ping-pong ..
	this.currentState = next;
}


SurfaceRenderer.prototype.resize = function(width, height)
{
	this.width = width;
	this.height = height;

	this.pathStates = [new SurfaceRendererState(this.width, this.height), 
					   new SurfaceRendererState(this.width, this.height)];

	this.reset();
}


