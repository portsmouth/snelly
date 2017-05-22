

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
	this.fbo == null;
	this.max_downres = 3;
	this.numSamples = 0;
	this.frametime_measure_ms = 0.0;
	this.spp = 0.0;

	// Default user-adjustable settings 
	this.renderMode = 'pt';
	this.maxBounces = 3;
	this.maxMarchSteps = 256;
	this.radianceClamp = 1.0e3;
	this.skyPower = 1.0;
	this.skyTemperature = 6000.0;
	this.exposure = 1.0;
	this.gamma = 2.2;
	this.whitepoint = 2.0;
	this.skipProbability = 0.0;
	this.goalFrametimeMs = 60.0;

	// Load shaders
	this.shaderSources = GLU.resolveShaderSource(["pathtracer", "tonemapper", "filter"]);
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
	//this.frametime_measure_ms = 0.0;

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
	if (typeof sceneObj.shader == "undefined") { GLU.fail('Scene must define a "shader" function!'); }
	var shader = sceneObj.shader();

	// Insert default functions for those not implemented 
	if (shader.indexOf("SDF_SURFACE(")                     == -1) { shader += `\n float SDF_SURFACE(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
 	if (shader.indexOf("SURFACE_DIFFUSE_REFLECTANCE(")     == -1) { shader += `\n vec3 SURFACE_DIFFUSE_REFLECTANCE(vec3 X) { return vec3(1.0); }\n`; }
 	if (shader.indexOf("SURFACE_SPECULAR_REFLECTANCE(")    == -1) { shader += `\n vec3 SURFACE_SPECULAR_REFLECTANCE(vec3 X) { return vec3(1.0); }\n`; }
 	if (shader.indexOf("SURFACE_IOR(")                     == -1) { shader += `\n float SURFACE_IOR(vec3 X) { return 1.0; }\n`; }
 	if (shader.indexOf("SURFACE_ROUGHNESS(")               == -1) { shader += `\n float SURFACE_ROUGHNESS(vec3 X) { return 1.0; }\n`; }
 	
 	if (shader.indexOf("SDF_METAL(")                       == -1) { shader += `\n float SDF_METAL(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
 	if (shader.indexOf("METAL_SPECULAR_REFLECTANCE(")      == -1) { shader += `\n vec3 METAL_SPECULAR_REFLECTANCE(vec3 X) { return vec3(1.0); }\n`; }
 	if (shader.indexOf("METAL_ROUGHNESS(")                 == -1) { shader += `\n float METAL_ROUGHNESS(vec3 X) { return 1.0; }\n`; }
 	
 	if (shader.indexOf("SDF_DIELECTRIC(")                  == -1) { shader += `\n float SDF_DIELECTRIC(vec3 X) { const float HUGE_VAL = 1.0e20; return HUGE_VAL; }\n`; }
 	if (shader.indexOf("DIELECTRIC_SPECULAR_REFLECTANCE(") == -1) { shader += `\n vec3 DIELECTRIC_SPECULAR_REFLECTANCE(vec3 X) { return vec3(1.0); }\n`; }
 	if (shader.indexOf("DIELECTRIC_ROUGHNESS(")            == -1) { shader += `\n float DIELECTRIC_ROUGHNESS(vec3 X) { return 1.0; }\n`; }

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
			replacements.ENTRY_AO = 'main';
			this.aoProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
			break;

		case 'normals':
			replacements.ENTRY_NORMALS = 'main';
			this.normalsProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
			break;

		case 'pt':
		default:
			replacements.ENTRY_PATHTRACE_ALL = 'main';
			replacements.ENTRY_PATHTRACE_BLOCKS = 'ENTRY_PATHTRACE_BLOCKS';
			this.pathtraceAllProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
			replacements.ENTRY_PATHTRACE_BLOCKS = 'main';
			replacements.ENTRY_PATHTRACE_ALL = 'ENTRY_PATHTRACE_ALL';
			this.pathtraceBlocksProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
			break;
	}

	// Compile filter programs for all required kernel widths (once only):
	if (this.renderMode == 'pt')
	{
		if (this.filterPrograms == null)
		{
			this.filterPrograms = {};
			for (var d=this.max_downres; d>=1; d--) 
			{
				replacements.DOWN_RES = d;
				this.filterPrograms[d] = new GLU.Shader('filter', this.shaderSources, replacements);
			}
		}
	}

	// Tonemapping program
	this.tonemapProgram = new GLU.Shader('tonemapper', this.shaderSources, null);
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
	if (!this.loaded) return;

	var sceneObj = snelly.getScene(); if (sceneObj==null) return;
	if (snelly.getSpectra()==null) return;

	var timer_start = performance.now();

	var gl = this.gl;
	gl.disable(gl.DEPTH_TEST);
	gl.viewport(0, 0, this.width, this.height);

	////////////////////////////////////////////////
	/// Pathtracing
	////////////////////////////////////////////////

	var INTEGRATOR_PROGRAM = null;
	var DOWN_RES = 1.0;

	switch (this.renderMode)
	{
		case 'ao':
			INTEGRATOR_PROGRAM = this.aoProgram;
			break;

		case 'normals':
			INTEGRATOR_PROGRAM = this.normalsProgram;
			break;

		case 'pt':
		default:
			// Choose pathtracing mode according to total sample count:
			var DOWN_RES = Math.floor(Math.sqrt(this.width * this.height / Math.max(this.numSamples, 1)));
			DOWN_RES = Math.max(1, Math.min(this.max_downres, DOWN_RES));
			if (this.frametime_measure_ms < 2.0*this.goalFrametimeMs) DOWN_RES = 1.0;
			if (DOWN_RES <= 1.0) 
			{
				// Pathtrace all pixels
				INTEGRATOR_PROGRAM = this.pathtraceAllProgram;
			}
			else
			{
				// Pathtrace a subset, until > 1 spp has accumulated
				INTEGRATOR_PROGRAM = this.pathtraceBlocksProgram;
			}
			break;
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
	INTEGRATOR_PROGRAM.uniformF("camNear", camera.near);
	INTEGRATOR_PROGRAM.uniformF("camFar", camera.far);
	INTEGRATOR_PROGRAM.uniformF("camFovy", camera.fov);
	INTEGRATOR_PROGRAM.uniformF("camZoom", camera.zoom);
	INTEGRATOR_PROGRAM.uniformF("camAspect", camera.aspect);
	INTEGRATOR_PROGRAM.uniform2Fv("resolution", [this.width, this.height]);

	var sceneScale = 1.0;
	if (typeof sceneObj.getScale !== "undefined") 
	{
		sceneScale = sceneObj.getScale();
	}
	INTEGRATOR_PROGRAM.uniformF("sceneScale", sceneScale); 

	// Read wavelength -> XYZ table
	snelly.wavelengthToXYZ.bind(2);
	INTEGRATOR_PROGRAM.uniformTexture("WavelengthToXYZ", snelly.wavelengthToXYZ);
	snelly.emissionIcdf.bind(3);
	INTEGRATOR_PROGRAM.uniformTexture("ICDF", snelly.emissionIcdf);
	
	// Pathtracing options
	INTEGRATOR_PROGRAM.uniformI("maxBounces", this.maxBounces);
	INTEGRATOR_PROGRAM.uniformI("downRes", DOWN_RES);
	INTEGRATOR_PROGRAM.uniformF("skyPower", this.skyPower);
	INTEGRATOR_PROGRAM.uniformF("radianceClamp", Math.pow(10.0, this.radianceClamp));
	INTEGRATOR_PROGRAM.uniformF("skipProbability", this.skipProbability);

	// gamma correction of env map already done if sRGB ext was loaded
	if (GLU.sRGBExt != null) INTEGRATOR_PROGRAM.uniformF("gamma", 1.0);
	else                     INTEGRATOR_PROGRAM.uniformF("gamma", this.gamma);

	this.fbo.bind();
	this.fbo.drawBuffers(2);
	var current = this.currentState;
	var next    = 1 - current;
	this.pathStates[current].bind(INTEGRATOR_PROGRAM); // Read data from the 'current' state
	this.pathStates[next].attach(this.fbo);            // Write data into the 'next' state

	INTEGRATOR_PROGRAM.uniformI("haveEnvMap", Boolean(this.envMap) ? 1 : 0);
	if (this.envMap != null)
	{
		gl.activeTexture(gl.TEXTURE0 + 7);
		gl.bindTexture(gl.TEXTURE_2D, this.envMap.tex);
		var id = gl.getUniformLocation(INTEGRATOR_PROGRAM.program, "envMap");
		gl.uniform1i(id, 7);
	}

	// Upload current scene shader parameters
	if (typeof sceneObj.syncShader !== "undefined") 
	{
		sceneObj.syncShader(INTEGRATOR_PROGRAM); 
	}
	snelly.materials.syncShader(INTEGRATOR_PROGRAM);

	// Trace one path per pixel
	gl.disable(gl.BLEND);
	this.quadVbo.bind();
	this.quadVbo.draw(INTEGRATOR_PROGRAM, gl.TRIANGLE_FAN);
	this.fbo.unbind();

	////////////////////////////////////////////////
	/// Filtering / progressive rendering
	////////////////////////////////////////////////

	if ((this.renderMode=='pt') && (DOWN_RES>1))
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
	if ((this.renderMode!='pt') || DOWN_RES <= 1) 
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

	// Frame timing
	var timer_end = performance.now();
	var frame_time_ms = (timer_end - timer_start);
	if (this.numSamples==0)
	{
		this.frametime_measure_ms = frame_time_ms;
	}
	else
	{
		var smoothing = 0.05;
		var prev_frame_time_ms = this.frametime_measure_ms;
		this.frametime_measure_ms = smoothing*frame_time_ms+ (1.0-smoothing)*prev_frame_time_ms;
	}

	// Update skip probability to try to achieve interactive framerate
	var goalMs = Math.max(1.0, this.goalFrametimeMs);
	this.skipProbability = Math.max(0.0, 1.0-goalMs/this.frametime_measure_ms);

	// Update sample count
	this.numSamples += (1.0-this.skipProbability) * this.width * this.height / (DOWN_RES * DOWN_RES);
	this.spp = this.numSamples / (this.width * this.height);

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
	this.reset(true);
}


