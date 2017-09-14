
var PathtracerState = function(width, height)
{
	var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count 
	var rngData      = new Float32Array(width*height*4); // Random number seed
	for (var i=0; i<width*height; ++i)
	{
		rngData[i*4 + 0] = Math.random()*4194167.0;
		rngData[i*4 + 1] = Math.random()*4194167.0;
		rngData[i*4 + 2] = Math.random()*4194167.0;
		rngData[i*4 + 3] = Math.random()*4194167.0;
	}
	this.radianceTex = new GLU.Texture(width, height, 4, true, false, true, radianceData);
	this.rngTex      = new GLU.Texture(width, height, 4, true, false, true, rngData);
}

PathtracerState.prototype.bind = function(shader)
{
	this.radianceTex.bind(0);
	this.rngTex.bind(1);
	shader.uniformTexture("Radiance", this.radianceTex);
	shader.uniformTexture("RngData", this.rngTex);
}

PathtracerState.prototype.attach = function(fbo)
{
	var gl = GLU.gl;
	fbo.attachTexture(this.radianceTex, 0);
	fbo.attachTexture(this.rngTex, 1);
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

/** 
* Interface to the renderer. The rendering modes available are:
*  - 'pt': pathtracer (uni-directional)
*  - 'ao': ambient occlusion, colored via {@link Surface} material diffuse albedo modulated by the `SURFACE_DIFFUSE_REFLECTANCE` shader function
*  - 'normals': view normal at first hit as a color
*  - 'firsthit': view first hit {@link Surface} material diffuse albedo, modulated by the sum of 
*    `SURFACE_DIFFUSE_REFLECTANCE` and `SURFACE_SPECULAR_REFLECTANCE` shader functions
* @constructor 
* @property {number} width                 - (if not specified, fits to window) 
* @property {number} height                - (if not specified, fits to window) 
* @property {String} [renderMode='pt']     - rendering mode (either 'pt', 'ao', 'normals', 'firsthit') 
* @property {number} [maxMarchSteps=256]   - maximum number of raymarching steps per path segment
* @property {number} [radianceClamp=3.0]   - clamp radiance to (10^) this max value, for firefly reduction
* @property {number} [skyPower=4.0]        - sky power (arbitrary units)
* @property {number} [skyTemperature=6000] - sky temperature (in Kelvin)
* @property {number} [exposure=4.5]        - exposure, on a log scale
* @property {number} [gamma=2.2]           - display gamma correction
* @property {number} [whitepoint=2.0]      - tonemapping whitepoint
* @property {number} [goalFPS=10.0]        - sampling will adjust to try to match goal FPS
* @property {number} [minsSPPToRedraw=0.0] - if >0.0, renderer will not redraw until the specified SPP have been accumulated
* @property {number} [envMapVisible=true]  - whether env map is visible to primary rays (otherwise black)
* @property {number} [envMapRotation=0.0]  - env map rotation about pole in degrees (0 to 360)
* @property {number} [shadowStrength=1.0]  - if <1.0, areas in shadow are not completely dark (provided mostly to allow rendering of occluded areas, e.g. fractals)
* @property {number} [maxStepsIsMiss=true] - whether rays which exceed max step count are considered hits or misses
* @property {number} [AA=true]             - whether to jitter primary ray within pixel for AA
*/
var Renderer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = snelly.render_canvas;
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this._width = render_canvas.width;
	this._height = render_canvas.height;

	// Initialize pathtracing buffers and programs
	this.currentState = 0;
	this.pathStates = [new PathtracerState(this._width, this._height), 
					   new PathtracerState(this._width, this._height)];
	this.fbo == null;
	this.aoProgram           = null;
	this.normalsProgram      = null;
	this.firsthitProgram     = null;
	this.pathtraceAllProgram = null;
	this.tonemapProgram      = null;
	this.pickProgram         = null;

	// Internal properties (@todo: use underscore to make this more explicit?)
	this.numSamples = 0;
	this.numFramesSinceReset = 0;
	this.numFramesSinceInit = 0;
	this.skipProbability = 0.0;
	this.frametime_measure_ms = 0.0;
	this.spp = 0.0;

	// Default user-adjustable properties
	this.renderMode = 'pt';
	this.maxBounces = 3;
	this.maxMarchSteps = 256;
	this.radianceClamp = 3.0;
	this.skyPower = 1.0;
	this.skyTemperature = 6000.0;
	this.exposure = 3.0;
	this.gamma = 2.2;
	this.whitepoint = 2.0;
	this.goalFPS = 20.0;
	this.minsSPPToRedraw = 0.0;
	this.envMapVisible = true;
	this.envMapRotation = 0.0;
	this.shadowStrength = 1.0;
	this.maxStepsIsMiss = true;
	this.AA = true;

	// Load shaders
	this.shaderSources = GLU.resolveShaderSource(["pathtracer", "ao", "normals", "firsthit", "pick", "tonemapper"]);
	this.filterPrograms = null;

	// load env map
	this.loaded = true;
	var sceneObj = snelly.getScene();
	this.envMap = null;
	if (typeof sceneObj.envMap !== "undefined")
  	{
  		var url = sceneObj.envMap();
  		if (url != "")
  		{
  			var pathtracer = this;
  			this.loaded = false;
  			(function() { GLU.loadImageAndCreateTextureInfo(url, 
					function(imgInfo)
					{
						pathtracer.loaded =  true;	
						pathtracer.envMap = imgInfo;
						//computeEnvMapCDF();
					});
  			})(pathtracer.loaded);
  		}
  	}

	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	this.quadVbo = this.createQuadVbo();
	this.fbo = new GLU.RenderTarget();
	
	// Trigger initial buffer generation
	this.resize(this._width, this._height);
}

Renderer.prototype.createQuadVbo = function()
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


Renderer.computeEnvMapCDF = function()
{
	if (this.envMap == null) return;



}

/**
* Restart accumulating samples.
* @param {Boolean} [no_recompile=false] - set to true if shaders need recompilation too
*/
Renderer.prototype.reset = function(no_recompile = false)
{
	this.numSamples = 0;
	this.numFramesSinceReset = 0;
	if (!no_recompile) this.compileShaders();
	this.currentState = 0;
	this.pathStates[this.currentState].clear(this.fbo);
	this.pathStates[this.currentState+1].clear(this.fbo);
}

/**
* Read access to the current per-pixel sample count average.
* @returns {number} - the current sample count average
*/
Renderer.prototype.getSPP = function()
{
	return this.spp;
}

Renderer.prototype.compileShaders = function()
{
	// Inject code for the current scene SDF:
	var sceneObj = snelly.getScene();
	if (sceneObj == null) return;
	if (typeof sceneObj.shader == "undefined") { GLU.fail('Scene must define a "shader" function!'); }
	var shader = sceneObj.shader();

	// Insert default functions for those not implemented 
	if (shader.indexOf("INIT(")                            == -1) { shader += `\n void INIT() {}\n`; }

	if (shader.indexOf("SDF_SURFACE(")                     == -1) { shader += `\n float SDF_SURFACE(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
 	if (shader.indexOf("SURFACE_DIFFUSE_REFLECTANCE(")     == -1) { shader += `\n vec3 SURFACE_DIFFUSE_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
 	if (shader.indexOf("SURFACE_SPECULAR_REFLECTANCE(")    == -1) { shader += `\n vec3 SURFACE_SPECULAR_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
 	if (shader.indexOf("SURFACE_ROUGHNESS(")               == -1) { shader += `\n float SURFACE_ROUGHNESS(in float roughness, in vec3 X, in vec3 N) { return roughness; }\n`; }
 	
 	if (shader.indexOf("SDF_METAL(")                       == -1) { shader += `\n float SDF_METAL(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
 	if (shader.indexOf("METAL_SPECULAR_REFLECTANCE(")      == -1) { shader += `\n vec3 METAL_SPECULAR_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
 	if (shader.indexOf("METAL_ROUGHNESS(")                 == -1) { shader += `\n float METAL_ROUGHNESS(in float roughness, in vec3 X, in vec3 N) { return roughness; }\n`; }
 	
 	if (shader.indexOf("SDF_DIELECTRIC(")                  == -1) { shader += `\n float SDF_DIELECTRIC(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
 	if (shader.indexOf("DIELECTRIC_SPECULAR_REFLECTANCE(") == -1) { shader += `\n vec3 DIELECTRIC_SPECULAR_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
 	if (shader.indexOf("DIELECTRIC_ROUGHNESS(")            == -1) { shader += `\n float DIELECTRIC_ROUGHNESS(in float roughness, in vec3 X, in vec3 N) { return roughness; }\n`; }

	if (shader.indexOf("SDF_VOLUME(")                      == -1) { shader += `\n float SDF_VOLUME(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
	if (shader.indexOf("VOLUME_EXTINCTION(")               == -1) { shader += `\n float VOLUME_EXTINCTION(vec3 X) { return 0.0; }\n`; }
	if (shader.indexOf("VOLUME_EXTINCTION_MAX(")           == -1) { shader += `\n float VOLUME_EXTINCTION_MAX() { return 0.0; }\n`; }
	if (shader.indexOf("VOLUME_ALBEDO(")                   == -1) { shader += `\n vec3 VOLUME_ALBEDO(vec3 X) { return vec3(0.0); }\n`; }
	if (shader.indexOf("VOLUME_EMISSION(")                 == -1) { shader += `\n vec3 VOLUME_EMISSION(vec3 X) { return vec3(0.0); }\n`; }

	// @todo: insert material GLSL code
	var dielectricObj = snelly.getLoadedDielectric(); if (dielectricObj == null) return;
	var metalObj      = snelly.getLoadedMetal();      if (metalObj == null) return;
	iorCodeDiele    = dielectricObj.ior();
	iorCodeMetal    = metalObj.ior();

	// Copy the current scene and material routines into the source code
	// of the trace fragment shader
	replacements = {};
	replacements.SHADER          = shader;
	replacements.IOR_FUNC        = iorCodeDiele + '\n' + iorCodeMetal;
	replacements.MAX_MARCH_STEPS = this.maxMarchSteps;
	replacements.MAX_BOUNCES     = this.maxBounces;

	// Compile pathtracer with different entry point according to mode.
	// Here shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	switch (this.renderMode)
	{
		case 'ao':
			this.aoProgram = new GLU.Shader('ao', this.shaderSources, replacements);
			break;
		case 'normals':
			this.normalsProgram = new GLU.Shader('normals', this.shaderSources, replacements);
			break;
		case 'firsthit':
			this.firsthitProgram = new GLU.Shader('firsthit', this.shaderSources, replacements);
			break;
		case 'pt':
		default:
			this.pathtraceAllProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
			break;
	}

	this.pickProgram = new GLU.Shader('pick', this.shaderSources, replacements);

	// Tonemapping program
	this.tonemapProgram = new GLU.Shader('tonemapper', this.shaderSources, null);
}

Renderer.prototype.enabled = function()
{
	return this.enable;
}

Renderer.prototype.depthTestEnabled = function()
{
	return this.depthTest;
}

Renderer.prototype.pick = function(xPick, yPick)
{
	if (!this.loaded) return;
	if (this.pickProgram==null) return;
	var sceneObj = snelly.getScene(); if (sceneObj==null) return;

	let gl = this.gl;
	gl.viewport(0, 0, 1, 1);

	let PROGRAM = this.pickProgram;
	PROGRAM.bind();

	// Upload current scene shader parameters
	if (typeof sceneObj.syncShader !== "undefined") 
	{
		sceneObj.syncShader(snelly, this.pickProgram); 
	}

	// sync camera info to shader	 
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camUp = camera.up.clone();
	camUp.transformDirection( camera.matrixWorld );
	var camX = new THREE.Vector3();
	camX.crossVectors(camUp, camDir);
	PROGRAM.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	PROGRAM.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	PROGRAM.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	PROGRAM.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
	PROGRAM.uniformF("camFovy", camera.fov);
	PROGRAM.uniformF("camAspect", camera.aspect);
	PROGRAM.uniformF("camAperture", camera.aperture);

	PROGRAM.uniformF("camFocalDistance", Math.max(1.0e-3*snelly.lengthScale, camera.focalDistance));
	PROGRAM.uniform2Fv("resolution", [this._width, this._height]);
	PROGRAM.uniformF("minLengthScale", snelly.minLengthScale);
	PROGRAM.uniformF("maxLengthScale", snelly.maxLengthScale);
	PROGRAM.uniformI("maxStepsIsMiss", Boolean(this.maxStepsIsMiss) ? 1 : 0);
	PROGRAM.uniform2Fv("mousePick", [xPick, yPick]); // picked pixel

	let fbo = new GLU.RenderTarget();
	fbo.bind();
	fbo.drawBuffers(1);

	let pickData = new Float32Array(4);
	this.pickTex = new GLU.Texture(1, 1, 4, true, false, true, pickData);
	fbo.attachTexture(this.pickTex, 0);

	// Trace pick ray
	this.quadVbo.bind();
	this.quadVbo.draw(PROGRAM, gl.TRIANGLE_FAN);

	// Read floating point hit distance output from pick fragment shader
	var pixels = new Float32Array(4);
	gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, pixels);
	fbo.unbind();

	return {distance: pixels[0], 
		    material: pixels[1]};
}

Renderer.prototype.render = function()
{
	if (!this.loaded) return;
	var sceneObj = snelly.getScene(); if (sceneObj==null) return;
	if (snelly.getSpectra()==null) return;

	var timer_start = performance.now();

	let gl = this.gl;
	if (typeof sceneObj.preframeCallback != "undefined")
	{
		sceneObj.preframeCallback(snelly, gl);
	}

	gl.disable(gl.DEPTH_TEST);
	gl.viewport(0, 0, this._width, this._height);

	// Update expected sample count after this frame, based on current resolution and skip probability
	this.numSamples += (1.0-this.skipProbability) * this._width * this._height;
	this.spp = this.numSamples / (this._width * this._height);

	// Update framebuffer only if we reached the requied SPP threshold for redraw
	if (this.spp >= this.minsSPPToRedraw)
	{
		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	}

	////////////////////////////////////////////////
	/// Pathtracing
	////////////////////////////////////////////////

	let INTEGRATOR_PROGRAM = null;
	switch (this.renderMode)
	{
		case 'ao':       INTEGRATOR_PROGRAM = this.aoProgram;           break;
		case 'normals':  INTEGRATOR_PROGRAM = this.normalsProgram;      break;
		case 'firsthit': INTEGRATOR_PROGRAM = this.firsthitProgram;     break;
		case 'pt':
		default:         INTEGRATOR_PROGRAM = this.pathtraceAllProgram; break;
	}

	INTEGRATOR_PROGRAM.bind();

	// sync camera info to shader	 
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camUp = camera.up.clone();
	camUp.transformDirection( camera.matrixWorld );
	var camX = new THREE.Vector3();
	camX.crossVectors(camUp, camDir);
	INTEGRATOR_PROGRAM.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	INTEGRATOR_PROGRAM.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	INTEGRATOR_PROGRAM.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	INTEGRATOR_PROGRAM.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
	INTEGRATOR_PROGRAM.uniformF("camFovy", camera.fov);
	INTEGRATOR_PROGRAM.uniformF("camAspect", camera.aspect);
	INTEGRATOR_PROGRAM.uniformF("camAperture", camera.aperture);
	INTEGRATOR_PROGRAM.uniformF("camFocalDistance", Math.max(1.0e-3*snelly.lengthScale, camera.focalDistance));
	INTEGRATOR_PROGRAM.uniform2Fv("resolution", [this._width, this._height]);

	// Read wavelength -> XYZ table
	if (this.renderMode=='pt' || this.renderMode=='ao' || this.renderMode=='firsthit')
	{
		snelly.wavelengthToXYZ.bind(2);
		INTEGRATOR_PROGRAM.uniformTexture("WavelengthToXYZ", snelly.wavelengthToXYZ);
		snelly.emissionIcdf.bind(3);
		INTEGRATOR_PROGRAM.uniformTexture("ICDF", snelly.emissionIcdf);
	}

	// Pathtracing options
	INTEGRATOR_PROGRAM.uniformI("maxBounces", this.maxBounces);
	INTEGRATOR_PROGRAM.uniformF("skyPower", this.skyPower);
	INTEGRATOR_PROGRAM.uniformF("radianceClamp", Math.pow(10.0, this.radianceClamp));
	INTEGRATOR_PROGRAM.uniformF("skipProbability", this.skipProbability);
	INTEGRATOR_PROGRAM.uniformF("maxLengthScale", snelly.maxLengthScale);
	INTEGRATOR_PROGRAM.uniformF("minLengthScale", snelly.minLengthScale);
	INTEGRATOR_PROGRAM.uniformF("shadowStrength", this.shadowStrength);
	INTEGRATOR_PROGRAM.uniformI("envMapVisible", Boolean(this.envMapVisible) ? 1 : 0);
	INTEGRATOR_PROGRAM.uniformF("envMapRotation", Math.min(Math.max(this.envMapRotation, 0.0), 360.0));
	INTEGRATOR_PROGRAM.uniformI("maxStepsIsMiss", Boolean(this.maxStepsIsMiss) ? 1 : 0);
	INTEGRATOR_PROGRAM.uniformF("jitter", Boolean(this.AA) ? 1 : 0);

	// Attach radiance FBO
	this.fbo.bind();
	this.fbo.drawBuffers(2);
	var current = this.currentState;
	var next    = 1 - current;
	this.pathStates[current].bind(INTEGRATOR_PROGRAM); // Read data from the 'current' state
	this.pathStates[next].attach(this.fbo);            // Write data into the 'next' state

	// Bind env map if we have it
	INTEGRATOR_PROGRAM.uniformI("haveEnvMap", Boolean(this.envMap) ? 1 : 0);
	if (this.envMap != null)
	{
		gl.activeTexture(gl.TEXTURE0 + 6);
		gl.bindTexture(gl.TEXTURE_2D, this.envMap.tex);
		var id = gl.getUniformLocation(INTEGRATOR_PROGRAM.program, "envMap");
		gl.uniform1i(id, 6);
	}

	// Upload current scene shader parameters
	if (typeof sceneObj.syncShader !== "undefined") 
	{
		sceneObj.syncShader(snelly, INTEGRATOR_PROGRAM); 
	}

	snelly.materials.syncShader(INTEGRATOR_PROGRAM);

	// Trace one path per pixel
	gl.disable(gl.BLEND);
	this.quadVbo.bind();
	this.quadVbo.draw(INTEGRATOR_PROGRAM, gl.TRIANGLE_FAN);
	this.fbo.unbind();

	////////////////////////////////////////////////
	/// Tonemapping / compositing
	////////////////////////////////////////////////

	if (this.spp >= this.minsSPPToRedraw)
	{
		this.tonemapProgram.bind();
		var radianceTexCurrent = this.pathStates[next].radianceTex;
		radianceTexCurrent.bind(0);	
		this.tonemapProgram.uniformTexture("Radiance", radianceTexCurrent);
		this.tonemapProgram.uniformF("exposure", this.exposure);
		this.tonemapProgram.uniformF("invGamma", 1.0/this.gamma);
		this.tonemapProgram.uniformF("whitepoint", this.whitepoint);

		gl.enable(gl.BLEND);
		gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
		gl.blendEquation(gl.FUNC_ADD);
		this.quadVbo.bind();
		this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);
	}

	// Ping-pong radiance buffers
	this.currentState = next;

	// Frame timing
	var timer_end = performance.now();
	var frame_time_ms = (timer_end - timer_start);
	if (this.numSamples==0)
	{
		this.frametime_measure_ms = frame_time_ms;
	}
	else
	{
		// apply exponential smoothing to frame time measurements
		var smoothing = 0.05;
		var prev_frame_time_ms = this.frametime_measure_ms;
		this.frametime_measure_ms = smoothing*frame_time_ms + (1.0-smoothing)*prev_frame_time_ms;
	}
	this.frametime_measure_ms = Math.max(1.0-6, this.frametime_measure_ms);

	// Update skip probability for next frame, to try to achieve interactive framerate
	let goalMs = Math.min(1.0e3, Math.max(1.0, 1.0e3/this.goalFPS));
	this.skipProbability = Math.max(0.0, 1.0-goalMs/this.frametime_measure_ms);

	this.numFramesSinceReset++;
	this.numFramesSinceInit++;

	if (typeof sceneObj.postframeCallback != "undefined")
	{
		sceneObj.postframeCallback(snelly, gl);
	}
}

Renderer.prototype.resize = function(width, height)
{
	this._width = width;
	this._height = height;

	this.fbo.unbind();
	// Two sets of radiance buffers, for read/write ping-pong
	this.pathStates = [new PathtracerState(this._width, this._height), 
					   new PathtracerState(this._width, this._height)];

	this.quadVbo = this.createQuadVbo();
	this.fbo = new GLU.RenderTarget();
	this.reset(true);
}

