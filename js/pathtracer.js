

var PathtracerState = function(width, height)
{
	var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count 
	var rngData      = new Float32Array(width*height*4); // Random number seed
	//var depthData    = new Float32Array(width*height*4); // Packed depth
	for (var i = 0; i<width*height; ++i)
	{
		for (var t = 0; t<4; ++t)
		{
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
	}
	this.radianceTex = new GLU.Texture(width, height, 4, true, false, true, radianceData);
	this.rngTex      = new GLU.Texture(width, height, 4, true, false, true, rngData);
	//this.depthTex    = new GLU.Texture(width, height, 4, true, false, true, depthData);
}

PathtracerState.prototype.bind = function(shader)
{
	this.radianceTex.bind(0);
	this.rngTex.bind(1);
	//this.depthTex.bind(2);
	shader.uniformTexture("Radiance", this.radianceTex);
	shader.uniformTexture("RngData", this.rngTex);
	//shader.uniformTexture("Depth", this.depthTex);
}

PathtracerState.prototype.attach = function(fbo)
{
	var gl = GLU.gl;
	fbo.attachTexture(this.radianceTex, 0);
	fbo.attachTexture(this.rngTex, 1);
	//fbo.attachTexture(this.depthTex, 2);
	if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) 
	{
		GLU.fail("Invalid framebuffer");
	}
}

PathtracerState.prototype.detach = function(fbo)
{
	var gl = GLU.gl;
	fbo.detachTexture(0);
	fbo.detachTexture(1);
	//fbo.detachTexture(2);
}

PathtracerState.prototype.clear = function(fbo)
{
	// clear radiance buffer
	var gl = GLU.gl;
	fbo.bind();
	fbo.drawBuffers(1);
	fbo.attachTexture(this.radianceTex, 0);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT);
	fbo.unbind();
}

var Pathtracer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;

	// Initialize pathtracing textures:
	this.currentState = 0;
	this.pathStates = [new PathtracerState(this.width, this.height), 
					   new PathtracerState(this.width, this.height)];

	this.maxBounces = 1;
	this.maxMarchSteps = 32;
	this.enable = true;
	//this.depthTest = false;
	this.showBounds = false;
	this.exposure = 10.0;
	this.fbo == null;
	this.max_downres = 4;
	this.numSamples = 0;
	this.skyPower = 1.0;
	this.diffuseAlbedo = [1.0, 1.0, 1.0];
	this.absorptionDiele = [0.0, 0.0, 0.0];

	// Load shaders
	this.shaderSources = GLU.resolveShaderSource(["pathtracer", "tonemapper", "pick", "filter"]);
	this.filterPrograms = null;
	this.compileShaders();

	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	this.quadVbo = this.createQuadVbo();
	this.fbo = new GLU.RenderTarget();
	
	// Trigger initial buffer generation
	this.resize(this.width, this.height);
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

Pathtracer.prototype.reset = function(no_recompile = false)
{
	this.downres = this.max_downres;
	this.numSamples = 0;

	if (!no_recompile) this.compileShaders();
	this.currentState = 0;
	this.pathStates[this.currentState].clear(this.fbo);
	this.pathStates[this.currentState+1].clear(this.fbo);
}


Pathtracer.prototype.compileShaders = function()
{
	// Inject code for the current scene SDF:
	var sceneObj = snelly.getScene();
	if (sceneObj == null) return;
	var sdfCode = sceneObj.sdf();

	// Insert dummy functions if missing 
	if (sdfCode.indexOf("SDF_METAL(")      == -1) { sdfCode += `\n float SDF_METAL(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
	if (sdfCode.indexOf("SDF_DIELECTRIC(") == -1) { sdfCode += `\n float SDF_DIELECTRIC(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
	if (sdfCode.indexOf("SDF_DIFFUSE(")    == -1) { sdfCode += `\n float SDF_DIFFUSE(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }

	// @todo: insert material GLSL code
	var dielectricObj = snelly.getLoadedDielectric(); if (dielectricObj == null) return;
	var metalObj      = snelly.getLoadedMetal();      if (metalObj == null) return;
	iorCodeDiele    = dielectricObj.ior();
	iorCodeMetal    = metalObj.ior();

	// Copy the current scene and material routines into the source code
	// of the trace fragment shader
	replacements = {};
	replacements.SDF_FUNC        = sdfCode;
	replacements.IOR_FUNC        = iorCodeDiele + '\n' + iorCodeMetal;
	replacements.MAX_MARCH_STEPS = this.maxMarchSteps;
	replacements.MAX_BOUNCES     = this.maxBounces;

	// Compile pathtracer with different entry point according to mode.
	// Here shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	replacements.RENDER_ALL = 'main';
	replacements.RENDER_BLOCKS = 'RENDER_BLOCKS';
	this.pathtraceAllProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);

	replacements.RENDER_BLOCKS = 'main';
	replacements.RENDER_ALL = 'RENDER_ALL';
	this.pathtraceBlocksProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);

	// Compile filter programs for all required kernel widths (once only):
	if (this.filterPrograms == null)
	{
		this.filterPrograms = {};
		for (var d=this.max_downres; d>=1; d--) 
		{
			replacements.DOWN_RES = d;
			this.filterPrograms[d] = new GLU.Shader('filter', this.shaderSources, replacements);
		}
	}

	// Adjunct programs
	this.pickProgram      = new GLU.Shader('pick',       this.shaderSources, replacements);
	this.tonemapProgram   = new GLU.Shader('tonemapper', this.shaderSources, null);
}


Pathtracer.prototype.pick = function(xPick, yPick)
{
	var sceneObj = snelly.getScene();
	if (sceneObj == null) return null;
	var gl = this.gl;
	gl.viewport(0, 0, 1, 1);

	this.pickProgram.bind();
	sceneObj.syncShader(this.pickProgram); 

	// sync camera info to shader	 
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camY = camera.up.clone();
	camY.transformDirection( camera.matrixWorld );
	var camX = new THREE.Vector3();
	camX.crossVectors(camY, camDir);

	var ndcX = xPick;
	var ndcY = yPick;

	this.pickProgram.uniformF("ndcX", ndcX);
	this.pickProgram.uniformF("ndcY", ndcY);
	this.pickProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	this.pickProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	this.pickProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	this.pickProgram.uniform3Fv("camY", [camY.x, camY.y, camY.z]);
	this.pickProgram.uniformF("camNear", camera.near);
	this.pickProgram.uniformF("camFar", camera.far);
	this.pickProgram.uniformF("camFovy", camera.fov);
	this.pickProgram.uniformF("camZoom", camera.zoom);
	this.pickProgram.uniformF("camAspect", camera.aspect);
	this.pickProgram.uniformF("SceneScale", sceneObj.getScale()); 

	var fbo = new GLU.RenderTarget();
	fbo.bind();
	fbo.drawBuffers(1);

	var pickData = new Float32Array(4);
	this.pickTex = new GLU.Texture(1, 1, 4, true, false, true, pickData);
	fbo.attachTexture(this.pickTex, 0);

	// Trace pick ray
	this.quadVbo.bind();
	this.quadVbo.draw(this.pickProgram, gl.TRIANGLE_FAN);

	// Read floating point hit distance output from pick fragment shader
	var pixels = new Uint8Array(4);
	gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
	pixels = new Float32Array(pixels.buffer);
	fbo.unbind();
	var hitDist = pixels[0];
	if (hitDist < 0.0) return null;

	// Compute hit point given pick NDC, camera, and the hit distance
	var fh = camera.near * Math.tan(0.5*camera.fov*Math.PI/180.0) / camera.zoom;
	var fw = camera.aspect * fh;

	var sX = camX.clone();
	var sY = camY.clone();
	sX.multiplyScalar(-fw*ndcX);
	sY.multiplyScalar(fh*ndcY);
	var s = sX.clone();
	s.add(sY);

	var D = camDir.clone();
	D.multiplyScalar(camera.near);
	D.add(s);
	D.normalize(); // ray direction through camPos, towards hit point
	D.multiplyScalar(hitDist);

	var hitPoint = camPos.clone();
	hitPoint.add(D);
	return hitPoint;
}

Pathtracer.prototype.enabled = function()
{
	return this.enable;
}

Pathtracer.prototype.depthTestEnabled = function()
{
	return this.depthTest;
}

Pathtracer.prototype.render = function()
{
	if (!this.enable) return;

	var sceneObj      = snelly.getScene();            if (sceneObj==null) return;
	var metalObj      = snelly.getLoadedMetal();      if (metalObj==null) return;
	var dielectricObj = snelly.getLoadedDielectric(); if (dielectricObj==null) return;
	if (snelly.getSpectra()==null) return;

	var gl = this.gl;
	gl.disable(gl.DEPTH_TEST);
	gl.viewport(0, 0, this.width, this.height);

	////////////////////////////////////////////////
	/// Pathtracing
	////////////////////////////////////////////////

	// Choose pathtracing mode according to total sample count:
	var DOWN_RES = Math.floor(Math.sqrt(this.width * this.height / Math.max(this.numSamples, 1)));
	DOWN_RES = Math.max(1, Math.min(this.max_downres, DOWN_RES));
	var PATHTRACER_PROGRAM = null;
	if (DOWN_RES <= 1.0) 
	{
		// Pathtrace all pixels
		PATHTRACER_PROGRAM = this.pathtraceAllProgram;
	}
	else
	{
		// Pathtrace a subset, until > 1 spp has accumulated
		PATHTRACER_PROGRAM = this.pathtraceBlocksProgram;
	}
	PATHTRACER_PROGRAM.bind();

	// sync camera info to shader	 
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camUp = camera.up.clone();
	camUp.transformDirection( camera.matrixWorld );
	var camX = new THREE.Vector3();
	camX.crossVectors(camUp, camDir);
	PATHTRACER_PROGRAM.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	PATHTRACER_PROGRAM.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	PATHTRACER_PROGRAM.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	PATHTRACER_PROGRAM.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
	PATHTRACER_PROGRAM.uniformF("camNear", camera.near);
	PATHTRACER_PROGRAM.uniformF("camFar", camera.far);
	PATHTRACER_PROGRAM.uniformF("camFovy", camera.fov);
	PATHTRACER_PROGRAM.uniformF("camZoom", camera.zoom);
	PATHTRACER_PROGRAM.uniformF("camAspect", camera.aspect);
	PATHTRACER_PROGRAM.uniform2Fv("resolution", [this.width, this.height]);
	PATHTRACER_PROGRAM.uniformF("SceneScale", sceneObj.getScale()); 

	// Read wavelength -> RGB table
	snelly.wavelengthToRgb.bind(2);
	PATHTRACER_PROGRAM.uniformTexture("WavelengthToRgb", snelly.wavelengthToRgb);
	snelly.emissionIcdf.bind(3);
	PATHTRACER_PROGRAM.uniformTexture("ICDF", snelly.emissionIcdf);
	
	// Pathtracing options
	PATHTRACER_PROGRAM.uniformI("MaxBounces", this.maxBounces);
	PATHTRACER_PROGRAM.uniformI("downRes", DOWN_RES);
	PATHTRACER_PROGRAM.uniformF("SkyPower", this.skyPower);
	PATHTRACER_PROGRAM.uniform3Fv("diffuseAlbedo", this.diffuseAlbedo);

	this.fbo.bind();
	this.fbo.drawBuffers(2);
	var current = this.currentState;
	var next    = 1 - current;
	this.pathStates[current].bind(PATHTRACER_PROGRAM); // Read data from the 'current' state
	this.pathStates[next].attach(this.fbo);            // Write data into the 'next' state

	// Upload current scene SDF shader parameters
	sceneObj.syncShader(PATHTRACER_PROGRAM); 
	dielectricObj.syncShader(PATHTRACER_PROGRAM); // upload current dielectric IOR parameters
	metalObj.syncShader(PATHTRACER_PROGRAM);      // upload current metal IOR parameters

	// Trace one path per pixel
	gl.disable(gl.BLEND);
	this.quadVbo.bind();
	this.quadVbo.draw(PATHTRACER_PROGRAM, gl.TRIANGLE_FAN);
	this.fbo.unbind();

	////////////////////////////////////////////////
	/// Filtering / progressive rendering
	////////////////////////////////////////////////

	if (DOWN_RES > 1) 
	{
		this.fbo.bind();
		this.fbo.drawBuffers(1);
		// Generate filtered image from sparse radiance samples
		this.filterPrograms[DOWN_RES].bind();
		this.pathStates[next].radianceTex.bind(0);    // read radiance
		this.filterPrograms[DOWN_RES].uniformTexture("Radiance", this.pathStates[next].radianceTex);
		this.fbo.attachTexture(this.filteredTex, 0);  // write filter-averaged mean
		this.filterPrograms[DOWN_RES].uniform2Fv("resolution", [this.width, this.height]);
		this.quadVbo.bind();
		this.quadVbo.draw(this.filterPrograms[DOWN_RES], gl.TRIANGLE_FAN);
		this.fbo.unbind();
	}

	////////////////////////////////////////////////
	/// Tonemapping / compositing
	////////////////////////////////////////////////

	this.tonemapProgram.bind();
	var radianceTexCurrent;
	if (DOWN_RES <= 1) 
	{
		// Use full resolution radiance buffer
		radianceTexCurrent = this.pathStates[next].radianceTex;
		//var depthTexCurrent = this.pathStates[next].depthTex;
	}
	else
	{
		// Use filtered radiance buffer until sufficient sample count reached
		radianceTexCurrent = this.filteredTex;
	}

	radianceTexCurrent.bind(0);	
	//depthTexCurrent.bind(1);
	this.tonemapProgram.uniformTexture("Radiance", radianceTexCurrent);
	//this.tonemapProgram.uniformTexture("DepthSurface", depthTexCurrent);
	this.tonemapProgram.uniformF("exposure", this.exposure);
	this.tonemapProgram.uniformF("invGamma", 1.0);
	this.tonemapProgram.uniformF("alpha", 1.0);

	gl.enable(gl.BLEND);
	gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
	gl.blendEquation(gl.FUNC_ADD);
	this.quadVbo.bind();
	this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);

	// Ping-pong ..
	this.currentState = next;

	// Update sample count
	this.numSamples += this.width * this.height / (DOWN_RES * DOWN_RES);

	this.fbo.unbind();
}


Pathtracer.prototype.resize = function(width, height)
{
	this.width = width;
	this.height = height;

	this.fbo.unbind();
	this.pathStates = [new PathtracerState(this.width, this.height), 
					   new PathtracerState(this.width, this.height)];

    // texture for progressive mode
	var filteredData = new Float32Array(width*height*4); 
	this.filteredTex = new GLU.Texture(width, height, 4, true, false, true, filteredData);

	this.quadVbo = this.createQuadVbo();
	this.fbo = new GLU.RenderTarget();
	this.reset();
}


