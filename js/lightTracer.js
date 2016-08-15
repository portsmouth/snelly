

var LightTracerState = function(size) 
{
	var posData   = new Float32Array(size*size*4); // ray position
	var dirData   = new Float32Array(size*size*4); // ray direction
	var rngData   = new Float32Array(size*size*4); // Random number seed
	var rgbData   = new Float32Array(size*size*4); // Ray color, and wavelength

	for (var i = 0; i<size*size; ++i)
	{
		for (var t = 0; t<4; ++t)
		{
			dirData[i*4 + t] = 0.0;
			rgbData[i*4 + t] = Math.random();
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
		dirData[i*4 + 0] = 1.0;
	}

	this.posTex   = new GLU.Texture(size, size, 4, true, false, true, posData);
	this.dirTex   = new GLU.Texture(size, size, 4, true, false, true, dirData);
	this.rngTex   = new GLU.Texture(size, size, 4, true, false, true, rngData);
	this.rgbTex   = new GLU.Texture(size, size, 4, true, false, true, rgbData);
}

LightTracerState.prototype.bind = function(shader)
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

LightTracerState.prototype.attach = function(fbo)
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

LightTracerState.prototype.detach = function(fbo)
{
	var gl = GLU.gl;
	fbo.detachTexture(0);
	fbo.detachTexture(1);
	fbo.detachTexture(2);
	fbo.detachTexture(3);
}


var LightTracer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	// Initialize textures containing ray states
	this.raySize = 128;
	this.maxMarchSteps = 512;
	this.maxPathLength = 32;
	this.enabled = true;
	this.initStates();

	// Create a quad VBO for rendering textures
	this.quadVbo = this.createQuadVbo();

	// Spectrum initialization
	this.spectra = {}
	this.SPECTRUM_SAMPLES = 256;
	this.spectrumObj = null;
	this.LAMBDA_MIN = 360.0;
    this.LAMBDA_MAX = 750.0;
	var wToRgb = wavelengthToRgbTable();
	this.wavelengthToRgb = new GLU.Texture(wToRgb.length/4, 1, 4, true,  true, true, wToRgb);
	this.emissionIcdf    = new GLU.Texture(4*this.SPECTRUM_SAMPLES, 1, 1, true, false, true, null);

	this.addSpectrum( new FlatSpectrum("flat", "Flat spectrum", 400.0, 700.0) );
	this.addSpectrum( new BlackbodySpectrum("blackbody", "Blackbody spectrum", 6000.0) );
	this.addSpectrum( new MonochromaticSpectrum("monochromatic", "Monochromatic spectrum", 650.0) ); 

	this.fbo == null;
	this.loadSpectrum("flat");

	// Initialize raytracing shaders
	this.shaderSources = GLU.resolveShaderSource(["init", "trace", "line", "linedepth", "comp", "pass"]);

	// Initialize GL
	this.fbo = new GLU.RenderTarget();
}


LightTracer.prototype.createQuadVbo = function()
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

LightTracer.prototype.reset = function()
{
	this.wavesTraced = 0;
	this.raysTraced = 0;
	this.pathsTraced = 0;
	this.pathLength = 0;

	this.compileShaders();

	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.fluenceBuffer, 0);
	this.gl.clearColor(0.0, 0.0, 0.0, 1.0);
	this.gl.clear(this.gl.COLOR_BUFFER_BIT);
	this.fbo.unbind();

	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.depthTex[0], 0);
	this.gl.clearColor(1.0, 0.0, 0.0, 1.0);
	this.gl.clear(this.gl.COLOR_BUFFER_BIT);
	this.fbo.unbind();

	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.depthTex[1], 0);
	this.gl.clearColor(1.0, 0.0, 0.0, 1.0);
	this.gl.clear(this.gl.COLOR_BUFFER_BIT);
	this.fbo.unbind();
}


LightTracer.prototype.resetActiveBlock = function()
{
	//this.activeBlock = 4;
	this.activeBlock = this.raySize;
}


LightTracer.prototype.compileShaders = function()
{
	var sceneObj = snelly.getLoadedScene();
	var materialObj = snelly.getLoadedMaterial();
	if (sceneObj==null || materialObj==null) return;
	
	// Inject code for the current scene and material:
	sdfCode    = sceneObj.sdf();
	iorCode    = materialObj.ior();
	sampleCode = materialObj.sample();

	// Copy the current scene and material routines into the source code
	// of the trace fragment shader
	replacements = {};
	replacements.SDF_FUNC    = sdfCode;
	replacements.IOR_FUNC    = iorCode;
	replacements.SAMPLE_FUNC = sampleCode;
	replacements.MAX_MARCH_STEPS = this.maxMarchSteps;

	// shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.traceProgram = new GLU.Shader('trace', this.shaderSources, replacements);
	this.initProgram  = new GLU.Shader('init',  this.shaderSources, null);
	this.lineProgram  = new GLU.Shader('line',  this.shaderSources, null);
	this.lineDepthProgram  = new GLU.Shader('linedepth',  this.shaderSources, null);
	this.compProgram  = new GLU.Shader('comp',  this.shaderSources, null);
	this.passProgram  = new GLU.Shader('pass',  this.shaderSources, null);
}


LightTracer.prototype.initStates = function()
{
	this.resetActiveBlock();
	this.raySize = Math.floor(this.raySize);
	this.rayCount = this.raySize*this.raySize;
	this.currentState = 0;
	this.rayStates = [new LightTracerState(this.raySize), new LightTracerState(this.raySize)];
		
	// Create the buffer of texture coordinates, which maps each drawn line
	// to its corresponding texture lookup.
	{
		var gl = GLU.gl;
		this.rayVbo = new GLU.VertexBuffer();
		this.rayVbo.addAttribute("TexCoord", 3, gl.FLOAT, false);
		this.rayVbo.init(this.rayCount*2);
		var vboData = new Float32Array(this.rayCount*2*3);
		for (var i=0; i<this.rayCount; ++i)
		{
			var u = ((i % this.raySize) + 0.5) / this.raySize;
			var v = (Math.floor(i/this.raySize) + 0.5) / this.raySize;
			vboData[i*6 + 0] = vboData[i*6 + 3] = u;
			vboData[i*6 + 1] = vboData[i*6 + 4] = v;
			vboData[i*6 + 2] = 0.0;
			vboData[i*6 + 5] = 1.0;
		}
		this.rayVbo.copy(vboData);
	}
}

// emission spectrum management
LightTracer.prototype.addSpectrum = function(spectrumObj)
{
	this.spectra[spectrumObj.getName()] = spectrumObj;
}

LightTracer.prototype.getSpectra = function()
{
	return this.spectra;
}

LightTracer.prototype.loadSpectrum = function(spectrumName)
{
	this.spectrumObj = this.spectra[spectrumName];
	var inverseCDF = this.spectrumObj.inverseCDF(this.LAMBDA_MIN, this.LAMBDA_MAX, this.SPECTRUM_SAMPLES);
	this.emissionIcdf.bind(0);
    this.emissionIcdf.copy(inverseCDF);
    if (this.fbo != null)
		this.reset();
}

LightTracer.prototype.getLoadedSpectrum = function()
{
	return this.spectrumObj;
}


LightTracer.prototype.composite = function()
{
	// Normalize and tonemap the fluence buffer to produce final pixels
	this.fluenceBuffer.bind(0);
	this.compProgram.bind();
	this.quadVbo.bind();
	this.compProgram.uniformTexture("Fluence", this.fluenceBuffer);

	// Normalize the emission by dividing by the total number of paths
	// (and also apply gamma correction)
	var gui = snelly.getGUI();
	this.compProgram.uniformF("invNumPaths", 1.0/Math.max(this.pathsTraced, 1));
	this.compProgram.uniformF("exposure", gui.lightTracerSettings.exposure);
	this.compProgram.uniformF("invGamma", 1.0/gui.lightTracerSettings.gamma);

	this.quadVbo.draw(this.compProgram, this.gl.TRIANGLE_FAN);
}


LightTracer.prototype.isEnabled = function()
{
	return this.enabled;
}

LightTracer.prototype.render = function()
{
	if (!this.enabled) return;

	var sceneObj = snelly.getLoadedScene();
	if (sceneObj==null) return;

	var materialObj = snelly.getLoadedMaterial();
	if (materialObj==null) return;

	var current = this.currentState;
	var next    = 1 - current;

	// trace light beams
	this.fbo.bind();
	var gl = GLU.gl;
	gl.viewport(0, 0, this.raySize, this.raySize);
	//gl.scissor(0, 0, this.raySize, this.activeBlock);
	//gl.enable(gl.SCISSOR_TEST);

	// We will write the next state's ray data
	this.fbo.drawBuffers(4);
	this.rayStates[next].attach(this.fbo);
	this.quadVbo.bind();

	// initialize emitter rays (the beginning of a 'wave')
	if (this.pathLength == 0)
	{
		this.initProgram.bind(); // Start all rays at emission point(s)
		this.rayStates[current].rngTex.bind(0); // Read random seed from the current state
		this.initProgram.uniformTexture("RngData", this.rayStates[current].rngTex);

		// Read wavelength -> RGB table 
		this.wavelengthToRgb.bind(1);
		this.emissionIcdf.bind(2);
		this.initProgram.uniformTexture("WavelengthToRgb", this.wavelengthToRgb);
		this.initProgram.uniformTexture("ICDF", this.emissionIcdf);
		
		// Emitter data
		var laser = snelly.laser;
		emitterPos = laser.getPoint();
		emitterDir = laser.getDirection();
		emitterRadius = laser.getEmissionRadius();
		emissionSpread = laser.getEmissionSpreadAngle();
		this.initProgram.uniform3F("EmitterPos", emitterPos.x, emitterPos.y, emitterPos.z);
		this.initProgram.uniform3F("EmitterDir", emitterDir.x, emitterDir.y, emitterDir.z);
		this.initProgram.uniformF("EmitterRadius", emitterRadius);
		this.initProgram.uniformF("EmitterSpread", emissionSpread);

		// Write emitted ray initial conditions into 'next' state
		this.quadVbo.draw(this.initProgram, gl.TRIANGLE_FAN);

		// Make this initial state be 'current'
		current = 1 - current;
		next    = 1 - next;

		// And we prepare to write into the 'next' state
		this.fbo.drawBuffers(4);
		this.rayStates[next].attach(this.fbo);
	}

	// Fire rays into scene, generating new ray pos/dir data in 'next' rayStates textures
	{
		this.traceProgram.bind();
		this.rayStates[current].bind(this.traceProgram);       // Use the current state as the initial conditions
		this.traceProgram.uniformF("SceneScale", sceneObj.getScale()); 
		sceneObj.syncShader(this.traceProgram);                // upload current scene SDF shader parameters
		materialObj.syncShader(this.traceProgram);             // upload current material IOR parameters
		this.quadVbo.draw(this.traceProgram, gl.TRIANGLE_FAN); // Generate the next ray state
		this.rayStates[next].detach(this.fbo);
	}

	// Read this data to draw the next 'wavefront' of bounces (i.e. line segments) into the wave buffer
	{
		gl.viewport(0, 0, this.width, this.height);
		gl.disable(gl.DEPTH_TEST);
		this.fbo.drawBuffers(1);
		this.fbo.attachTexture(this.waveBuffer, 0); // write to wave buffer
		if (this.pathLength == 0 || this.wavesTraced==0)
		{
			// Clear wavebuffer before the first bounce
			gl.clearColor(0.0, 0.0, 0.0, 1.0);
			gl.clear(gl.COLOR_BUFFER_BIT);
		}

		// The float radiance channels of all lines are simply added each pass
		gl.blendFunc(gl.ONE, gl.ONE); // accumulate line segment radiances
		gl.enable(gl.BLEND);
		this.lineProgram.bind();

		// Setup projection matrix
		var camera = snelly.camera;
		var projectionMatrix = camera.projectionMatrix.toArray();
		var projectionMatrixLocation = this.lineProgram.getUniformLocation("u_projectionMatrix");
		gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

		// Setup modelview matrix (to match camera)
		camera.updateMatrixWorld();
		var matrixWorldInverse = new THREE.Matrix4();
		matrixWorldInverse.getInverse( camera.matrixWorld );
		var modelViewMatrix = matrixWorldInverse.toArray();
		var modelViewMatrixLocation = this.lineProgram.getUniformLocation("u_modelViewMatrix");
		gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

		this.rayStates[current].posTex.bind(0); // read PosDataA = current.posTex
		this.rayStates[   next].posTex.bind(1); // read PosDataB = next.posTex
		this.rayStates[current].rgbTex.bind(2); // read current  = current.rgbTex
		this.lineProgram.uniformTexture("PosDataA", this.rayStates[current].posTex);
		this.lineProgram.uniformTexture("PosDataB", this.rayStates[   next].posTex);
		this.lineProgram.uniformTexture("RgbData",  this.rayStates[current].rgbTex);

		this.rayVbo.bind(); // Binds the TexCoord attribute
		this.rayVbo.draw(this.lineProgram, gl.LINES, this.raySize*this.activeBlock*2);

		this.raysTraced += this.raySize*this.activeBlock;
		this.pathLength += 1;
		//this.fbo.detachTexture(0);

		// Line depth pass
		{
			gl.disable(gl.BLEND);
			//gl.blendFunc(gl.ONE, gl.ZERO); // overwrite depth buffer each bounce

			this.lineDepthProgram.bind();

			// Setup projection matrix
			var camera = snelly.camera;
			var projectionMatrix = camera.projectionMatrix.toArray();
			var projectionMatrixLocation = this.lineDepthProgram.getUniformLocation("u_projectionMatrix");
			gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

			// Setup modelview matrix (to match camera)
			camera.updateMatrixWorld();
			var matrixWorldInverse = new THREE.Matrix4();
			matrixWorldInverse.getInverse( camera.matrixWorld );
			var modelViewMatrix = matrixWorldInverse.toArray();
			var modelViewMatrixLocation = this.lineDepthProgram.getUniformLocation("u_modelViewMatrix");
			gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

			this.rayStates[current].posTex.bind(0); // read PosDataA = current.posTex
			this.rayStates[   next].posTex.bind(1); // read PosDataB = next.posTex
			this.depthTex[current].bind(2);         // read Depth = depthTex[current]

			this.lineDepthProgram.uniformTexture("PosDataA", this.rayStates[current].posTex);
			this.lineDepthProgram.uniformTexture("PosDataB", this.rayStates[   next].posTex);
			this.lineDepthProgram.uniformTexture("Depth",  this.depthTex[current]);
			this.lineDepthProgram.uniformF("camNear", camera.near);
			this.lineDepthProgram.uniformF("camFar", camera.far);
			
			this.fbo.drawBuffers(1);
			this.fbo.attachTexture(this.depthTex[next], 0); // write to depthTex[next]

			this.rayVbo.bind(); // Binds the TexCoord attribute
			this.rayVbo.draw(this.lineDepthProgram, gl.LINES, this.raySize*this.activeBlock*2);

			//this.fbo.detachTexture(0);
		}
	}

	if (this.pathLength==this.maxPathLength || this.wavesTraced<this.maxPathLength)
	{
		// Add the wavebuffer contents, a complete set of rendered path segments,
		// into the fluence buffer
		gl.enable(gl.BLEND);
		gl.blendFunc(gl.ONE, gl.ONE); // accumulate radiances of all line segments from current 'wave' of bounces
		this.fbo.attachTexture(this.fluenceBuffer, 0);
		this.waveBuffer.bind(0);
		this.passProgram.bind();
		this.passProgram.uniformTexture("Frame", this.waveBuffer);
		this.quadVbo.bind();
		this.quadVbo.draw(this.passProgram, gl.TRIANGLE_FAN);
		this.pathsTraced += this.raySize*this.activeBlock;
		if (this.pathLength == this.maxPathLength)
		{
			this.wavesTraced += 1;
			this.pathLength = 0;
		}
	}

	//this.activeBlock = Math.min(512, this.activeBlock + 4);

	this.fbo.unbind();
	gl.disable(gl.BLEND);

	// Final composite of normalized fluenceBuffer to window
	this.composite();

	// Update raytracing state
	this.currentState = next;
}


LightTracer.prototype.getDepthTexture = function()
{
	return this.depthTex[this.currentState];
}


LightTracer.prototype.resize = function(width, height)
{
	this.width = width;
	this.height = height;

	this.fluenceBuffer = new GLU.Texture(width, height, 4, true, false, true, null);
	this.waveBuffer    = new GLU.Texture(width, height, 4, true, false, true, null);
	this.depthTex      = [new GLU.Texture(width, height, 4, true, false, true, null), 
                          new GLU.Texture(width, height, 4, true, false, true, null)];
	this.resetActiveBlock();
	this.reset();
}

