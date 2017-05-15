

var PathtracerState = function(width, height)
{
	var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count 
	var rngData      = new Float32Array(width*height*4); // Random number seed
	for (var i = 0; i<width*height; ++i)
	{
		for (var t = 0; t<4; ++t)
		{
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
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

	// @todo: these setting are defaults which should be specifiable in the scene
	this.maxBounces = 3;
	this.maxMarchSteps = 128;
	this.enable = true;
	this.exposure = 10.0;
	this.gamma = 2.2;
	this.whitepoint = 2.0;
	this.fbo == null;
	this.max_downres = 1;
	this.numSamples = 0;
	this.skyPower = 1.0;
	
	this.surfaceDiffuseAlbedoRGB = [1.0, 1.0, 1.0];
	this.surfaceDiffuseAlbedoXYZ = rgbToXyz(this.surfaceDiffuseAlbedoRGB);
	this.surfaceSpecAlbedoRGB = [1.0, 1.0, 1.0];
	this.surfaceSpecAlbedoXYZ = rgbToXyz(this.surfaceSpecAlbedoRGB);
	this.surfaceRoughness = 0.1;
	this.surfaceIor = 1.5;
	this.absorptionDieleRGB = [0.0, 0.0, 0.0];

	// Load shaders
	this.shaderSources = GLU.resolveShaderSource(["pathtracer", "tonemapper", "filter"]);
	this.filterPrograms = null;
	this.compileShaders();

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
					});
  			})(pathtracer.loaded);
  		}
  	}

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
	if (sdfCode.indexOf("SDF_SURFACE(")    == -1) { sdfCode += `\n float SDF_SURFACE(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }

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
	this.tonemapProgram   = new GLU.Shader('tonemapper', this.shaderSources, null);
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
	if (!this.loaded) return;

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

	// Read wavelength -> XYZ table
	snelly.wavelengthToXYZ.bind(2);
	PATHTRACER_PROGRAM.uniformTexture("WavelengthToXYZ", snelly.wavelengthToXYZ);
	snelly.emissionIcdf.bind(3);
	PATHTRACER_PROGRAM.uniformTexture("ICDF", snelly.emissionIcdf);
	
	// Pathtracing options
	PATHTRACER_PROGRAM.uniformI("MaxBounces", this.maxBounces);
	PATHTRACER_PROGRAM.uniformI("downRes", DOWN_RES);
	PATHTRACER_PROGRAM.uniformF("SkyPower", this.skyPower);

	// gamma correction of env map already done if sRGB ext was loaded
	if (GLU.sRGBExt != null) PATHTRACER_PROGRAM.uniformF("gamma", 1.0);
	else                     PATHTRACER_PROGRAM.uniformF("gamma", this.gamma);

	// Material parameters (or multipliers on user defined values)
	PATHTRACER_PROGRAM.uniform3Fv("surfaceDiffuseAlbedoXYZ", this.surfaceDiffuseAlbedoXYZ);
	PATHTRACER_PROGRAM.uniform3Fv("surfaceSpecAlbedoXYZ", this.surfaceSpecAlbedoXYZ);
	PATHTRACER_PROGRAM.uniformF("surfaceRoughness", this.surfaceRoughness);
	PATHTRACER_PROGRAM.uniformF("surfaceIor", this.surfaceIor);

	this.fbo.bind();
	this.fbo.drawBuffers(2);
	var current = this.currentState;
	var next    = 1 - current;
	this.pathStates[current].bind(PATHTRACER_PROGRAM); // Read data from the 'current' state
	this.pathStates[next].attach(this.fbo);            // Write data into the 'next' state

	PATHTRACER_PROGRAM.uniformI("haveEnvMap", Boolean(this.envMap) ? 1 : 0);
	if (this.envMap != null)
	{
		gl.activeTexture(gl.TEXTURE0 + 7);
		gl.bindTexture(gl.TEXTURE_2D, this.envMap.tex);
		var id = gl.getUniformLocation(PATHTRACER_PROGRAM.program, "envMap");
		gl.uniform1i(id, 7);
	}


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
	}
	else
	{
		// Use filtered radiance buffer until sufficient sample count reached
		radianceTexCurrent = this.filteredTex;
	}

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


