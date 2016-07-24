

//////////////////////////////////////////////////////////////////////
// RayState
//////////////////////////////////////////////////////////////////////

var RayState = function(size) 
{
	this.size = size;

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


//////////////////////////////////////////////////////////////////////
// Renderer
//////////////////////////////////////////////////////////////////////

var Renderer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;
	
	this.scenes = {}
	this.materials = {}

	// Initialize THREE.js camera
	var VIEW_ANGLE = 45;
	var ASPECT = this.width / this.height ;
	var NEAR = 0.05;
	var FAR = 1000;
	this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
	this.camera.position.z = 50;
	this.camera.updateProjectionMatrix();
	
	// Setup three.js GL viewport renderer
	var ui_canvas = document.getElementById('ui-canvas');
	ui_canvas.style.top = 0;
	ui_canvas.style.position = 'fixed' 

	this.glRenderer = new THREE.WebGLRenderer( { canvas: ui_canvas,
											     alpha: true,
											     antialias: true } );
	this.glRenderer.setClearColor( 0x000000, 0 ); // the default
	this.glRenderer.setSize(this.width, this.height);
	this.glScene = new THREE.Scene();
	this.glScene.add(this.camera);

	var pointLight = new THREE.PointLight(0xa0a0a0);
	pointLight.position.x = 10;
	pointLight.position.y = 50;
	pointLight.position.z = 130;
	this.glScene.add(pointLight);

	var light = new THREE.AmbientLight( 0x808080 ); // soft white light
	this.glScene.add( light );

	this.container = document.getElementById('container');
	{
		this.stats = new Stats();
		this.stats.domElement.style.position = 'absolute';
		this.stats.domElement.style.top = '0px';
		this.container.appendChild( this.stats.domElement );
	}

	// Create user control system for camera
	this.controls = new THREE.OrbitControls(this.camera, this.glRenderer.domElement);
	this.controls.zoomSpeed = 0.5;
	this.controls.addEventListener( 'change', camChanged );

	// Setup Laser pointer
	this.laser = new LaserPointer(this.glRenderer, this.glScene, this.camera, this.controls);
	this.laser.setPosition(new THREE.Vector3(-5.0, 0.0, 0.0));
	this.laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));

	// Setup keypress events
	renderer = this;
	window.addEventListener('keydown', function(event) 
	{
		var charCode = (event.which) ? event.which : event.keyCode;
		switch (charCode)
		{
			case 122: // F11 key: go fullscreen
				var element	= document.body;
				if      ( 'webkitCancelFullScreen' in document ) element.webkitRequestFullScreen();
				else if ( 'mozCancelFullScreen'    in document ) element.mozRequestFullScreen();
				else console.assert(false);
				break;
			case 70: // F key: focus on emitter
				renderer.controls.object.zoom = renderer.controls.zoom0;
				renderer.controls.target.copy(renderer.laser.getPoint());
				renderer.controls.update();
				break; 
		}
	}, false);
	
	this.glRenderer.domElement.addEventListener( 'mousemove', this, false );
	this.glRenderer.domElement.addEventListener( 'mousedown', this, false );
	this.glRenderer.domElement.addEventListener( 'mouseup',   this, false );


	////////////////////////////////////////////////////////////
	// Initialize ray rendering
	////////////////////////////////////////////////////////////

	// Initialize textures containing ray states
	this.raySize = 64;
	this.maxPathLength = 20;
	this.initRayStates();

	// Create a quad VBO for rendering textures
	this.quadVbo = this.createQuadVbo();

	// Spectrum initialization:
	// table of 256 vec4 RGB colors, corresponding to the 256 wavelength samples
	// between 360.0 and 750.0 nanometres
	// @todo:  for now we will just assume a flat emission spectrum, i.e pure white light.
	this.LAMBDA_MIN = 360.0;
    this.LAMBDA_MAX = 750.0;
	this.spectrumTable = wavelengthToRgbTable();
	this.spectrum = new GLU.Texture(this.spectrumTable.length/4, 1, 4, true,  true, true, this.spectrumTable);

	// Initialize raytracing shaders
	this.shaderSources = GLU.resolveShaderSource(["init", "trace", "line", "comp", "pass"]);

	// Initialize GL
	this.fbo = new GLU.RenderTarget();
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);
	this.resize(this.width, this.height);
}


Renderer.prototype.handleEvent = function(event)
{
	switch (event.type)
	{
		case 'mousemove': this.onDocumentMouseMove(event); break;
		case 'mousedown': this.onDocumentMouseDown(event); break;
		case 'mouseup':   this.onDocumentMouseUp(event);   break;
	}
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

Renderer.prototype.reset = function()
{
	if (!this.needsReset) return;

	this.needsReset = false;
	this.wavesTraced = 0;
	this.raysTraced = 0;
	this.samplesTraced = 0;
	this.pathLength = 0;

	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.screenBuffer, 0);
	this.gl.clear(this.gl.COLOR_BUFFER_BIT);
	this.fbo.unbind();
}


Renderer.prototype.resetActiveBlock = function()
{
	//this.activeBlock = 4;
	this.activeBlock = this.raySize;
}

Renderer.prototype.compileShaders = function()
{
	// @todo:  here, paste the current sdf and sample routine into the shader source:
	sdfCode    = this.sceneObj.sdf();
	iorCode    = this.materialObj.ior();
	sampleCode = this.materialObj.sample();
	injectedCode = sdfCode + '\n' + iorCode + '\n' + sampleCode;

	console.log('code to inject: ' + injectedCode);

	// shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.initProgram  = new GLU.Shader('init',  this.shaderSources);
	this.traceProgram = new GLU.Shader('trace', this.shaderSources);
	this.lineProgram  = new GLU.Shader('line',  this.shaderSources);
	this.compProgram  = new GLU.Shader('comp',  this.shaderSources);
	this.passProgram  = new GLU.Shader('pass',  this.shaderSources);
}


Renderer.prototype.initRayStates = function()
{
	this.resetActiveBlock();
	this.rayCount = this.raySize*this.raySize;
	this.currentState = 0;
	this.needsReset = true;
	
	this.rayStates = [new RayState(this.raySize), new RayState(this.raySize)];
		
	// Create the buffer of texture coordinates, which maps each drawn line
	// to its corresponding texture lookup.
	{
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


// scene management
Renderer.prototype.addScene = function(sceneObj)
{
	this.scenes[sceneObj.getName()] = sceneObj;
}

Renderer.prototype.getScenes = function()
{
	return this.scenes;
}

Renderer.prototype.loadScene = function(sceneName)
{
	this.sceneObj = this.scenes[sceneName];
	this.compileShaders();
	this.sceneObj.setCam(this.controls, this.camera);
	this.sceneObj.setLaser(this.laser);
}


// material management
Renderer.prototype.addMaterial = function(materialObj)
{
	this.materials[materialObj.getName()] = materialObj;
}

Renderer.prototype.getMaterials = function()
{
	return this.materials;
}

Renderer.prototype.loadMaterial = function(materialName)
{
	this.materialObj = this.materials[materialName];
	this.compileShaders();
}



Renderer.prototype.composite = function()
{
	this.screenBuffer.bind(0);
	this.compProgram.bind();
	this.quadVbo.bind();
	this.compProgram.uniformTexture("Frame", this.screenBuffer);

	// Tonemap to effectively divide by total number of emitted photons
	// (and also apply gamma correction)
	var numPhotons = Math.max(this.samplesTraced, 1);
	this.compProgram.uniformF("Exposure", 0.01 * renderer.params.renderSettings.exposure / numPhotons);
	
	this.quadVbo.draw(this.compProgram, this.gl.TRIANGLE_FAN);
}


Renderer.prototype.render = function()
{
	this.controls.update();

	var gl = this.gl;
	this.needsReset = true;

	var current = this.currentState;
	var next    = 1 - current;

	// Render laser pointer
	gl.viewport(0, 0, this.width, this.height);
	this.laser.render();
	
	////////////////////////////////
	// trace light beams
	////////////////////////////////

	this.fbo.bind();
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
		this.spectrum.bind(1);
		this.initProgram.uniformTexture("Spectrum", this.spectrum);
		
		// Emitter data
		emitterPos = this.laser.getPoint();
		emitterDir = this.laser.getDirection();
		emitterRadius = this.laser.getEmissionRadius();
		emissionSpread = this.laser.getEmissionSpreadAngle();
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

	// Raytrace into scene, generating new ray pos/dir data in 'next' rayStates textures
	{
		this.traceProgram.bind();
		this.rayStates[current].bind(this.traceProgram);       // Use the current state as the initial conditions
		this.sceneObj.syncShader(this.traceProgram);           // set current scene parameters 
		this.materialObj.syncShader(this.traceProgram);        // set current material parameters 
		this.quadVbo.draw(this.traceProgram, gl.TRIANGLE_FAN); // Generate the next ray state
		this.rayStates[next].detach(this.fbo);
	}

	// Draw the next set of lines into the wave buffer
	{
		//gl.disable(gl.SCISSOR_TEST);
		gl.viewport(0, 0, this.width, this.height);
		this.fbo.drawBuffers(1);
		this.fbo.attachTexture(this.waveBuffer, 0);
		if (this.pathLength == 0 || this.wavesTraced==0)
		{
			// Clear wavebuffer before the first bounce
			gl.clear(gl.COLOR_BUFFER_BIT);
		}

		gl.enable(gl.BLEND);
		this.lineProgram.bind();

		// Setup projection matrix
		var projectionMatrix = this.camera.projectionMatrix.toArray();
		var projectionMatrixLocation = this.lineProgram.getUniformLocation("u_projectionMatrix");
		gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

		// Setup modelview matrix (to match camera)
		this.camera.updateMatrixWorld();
		var matrixWorldInverse = new THREE.Matrix4();
		matrixWorldInverse.getInverse( this.camera.matrixWorld );
		var modelViewMatrix = matrixWorldInverse.toArray();
		var modelViewMatrixLocation = this.lineProgram.getUniformLocation("u_modelViewMatrix");
		gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

		this.rayStates[current].posTex.bind(0); // PosDataA = current.posTex
		this.rayStates[   next].posTex.bind(1); // PosDataB = next.posTex
		this.rayStates[current].rgbTex.bind(2); // current  = current.rgbTex

		this.lineProgram.uniformTexture("PosDataA", this.rayStates[current].posTex);
		this.lineProgram.uniformTexture("PosDataB", this.rayStates[   next].posTex);
		this.lineProgram.uniformTexture("RgbData",  this.rayStates[current].rgbTex);

		this.rayVbo.bind(); // Binds the TexCoord attribute
		this.rayVbo.draw(this.lineProgram, gl.LINES, this.raySize*this.activeBlock*2);

		this.raysTraced += this.raySize*this.activeBlock;
		this.pathLength += 1;
	}

	this.quadVbo.bind();

	// Update the screenBuffer with the waveBuffer contents
	if (this.pathLength==this.maxPathLength || this.wavesTraced<this.maxPathLength)
	{
		this.fbo.attachTexture(this.screenBuffer, 0);
		this.waveBuffer.bind(0);
		this.passProgram.bind();
		this.passProgram.uniformTexture("Frame", this.waveBuffer);
		this.quadVbo.draw(this.passProgram, gl.TRIANGLE_FAN);
		this.samplesTraced += this.raySize*this.activeBlock;
		if (this.pathLength == this.maxPathLength)
		{
			this.wavesTraced += 1;
			this.pathLength = 0;
		}
	}

	gl.disable(gl.BLEND);
	this.fbo.unbind();

	// Final composite of screenBuffer to window
	this.composite();

	// Update raytracing state
	this.currentState = next;

	this.stats.update();
}


Renderer.prototype.onDocumentMouseMove = function(event)
{
	event.preventDefault();
	if (this.laser.onMouseMove(event)) this.reset();
}

Renderer.prototype.onDocumentMouseDown = function(event)
{
	event.preventDefault();
	this.laser.onMouseDown(event);
}

Renderer.prototype.onDocumentMouseUp = function(event)
{
	event.preventDefault();
	this.laser.onMouseUp(event);
}

Renderer.prototype.getCamera = function()
{
	return this.camera;
}

Renderer.prototype.resize = function(width, height)
{
	this.width = width;
	this.height = height;
	this.camera.aspect = width / height;
	this.camera.updateProjectionMatrix();
	this.screenBuffer = new GLU.Texture(this.width, this.height, 4, true, false, true, null);
	this.waveBuffer   = new GLU.Texture(this.width, this.height, 4, true, false, true, null);
	this.resetActiveBlock();
	this.reset();

	this.glRenderer.setSize(width, height);
}



