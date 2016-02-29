


//////////////////////////////////////////////////////////////////////
// Laser pointer UI
//////////////////////////////////////////////////////////////////////

var LaserPointer = function(glRenderer, glScene, glCamera) 
{
	var group = new THREE.Object3D();

	var housingGeo      = new THREE.CylinderGeometry( 0.2, 0.2, 2, 32 );
	var housingMaterial = new THREE.MeshPhongMaterial( { color: 0xdddddd, specular: 0x999999, shininess: 60 } );
	var housingObj = new THREE.Mesh( housingGeo, housingMaterial  );
	group.add(housingObj);

	var azimuthHandleGeo      = new THREE.TorusGeometry( 0.5, 0.12, 32, 100 );
	var azimuthHandleMaterial = new THREE.MeshPhongMaterial( { color: 0x000099, specular: 0x000099, shininess: 30 } );
	var azimuthHandleObj = new THREE.Mesh( azimuthHandleGeo, azimuthHandleMaterial  );
	//azimuthHandleObj.rotation.y = 0.5*Math.PI;
	group.add(azimuthHandleObj);

	var angleHandleGeo      = new THREE.TorusGeometry( 0.5, 0.12, 32, 100 );
	var angleHandleMaterial = new THREE.MeshPhongMaterial( { color: 0x999900, specular: 0x999900, shininess: 30 } );
	var angleHandleObj = new THREE.Mesh( angleHandleGeo, angleHandleMaterial  );
	azimuthHandleObj.rotation.y = 0.5*Math.PI;
	group.add(angleHandleObj);

	// proxy to represent the emission point of the laser
	var proxyGeo = new THREE.SphereGeometry(0.05);
	var proxyMaterial = new THREE.MeshPhongMaterial( { color: 0xff0000, specular: 0x999900, shininess: 10 } );
	var proxyObj = new THREE.Mesh( proxyGeo, proxyMaterial  );
	proxyObj.position.y = 1.0;
	group.add(proxyObj);
	group.quaternion.setFromAxisAngle( new THREE.Vector3( 0, 0, 1 ), -Math.PI / 2 );
	group.position.x = -10.0;

	glScene.add(group);

	// back plane to implement dragging 
	var backPlaneObj = new THREE.Mesh(
					new THREE.PlaneBufferGeometry( 2000, 2000, 8, 8 ),
					new THREE.MeshBasicMaterial( { visible: false } )
				);
	glScene.add(backPlaneObj);

	this.objects = {
		"housing": housingObj,
		"azimuthHandle": azimuthHandleObj,
		"angleHandle": angleHandleObj,
		"proxy": proxyObj,
		"backPlane": backPlaneObj,
		"group": group
	};

	this.glScene = glScene;
	this.glRenderer = glRenderer;
	this.glCamera = glCamera;
}

LaserPointer.prototype.getObjects = function()
{
	return this.objects;
} 

LaserPointer.prototype.render = function()
{
	this.glRenderer.render(this.glScene, this.glCamera);
}

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
	
	this.initialized = false;

	// Initialize THREE.js camera
	{
		var VIEW_ANGLE = 45;
		var ASPECT = this.width / this.height ;
		var NEAR = 0.1;
		var FAR = 10000;
		this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
		this.camera.position.z = 50;
		this.camera.updateProjectionMatrix();
	}

	var renderer = this;
	GLU.resolveShaderSource(["init", "trace", "line", "comp", "pass"],
		function onShadersLoaded(shaderSources)
		{
			if (!renderer.initialized)
			{
				renderer.setup(shaderSources);
				renderer.initialized = true;
			}
		}
	);

	////////////////////////////////////////////////////////////
	// Setup three.js GL renderer
	////////////////////////////////////////////////////////////

	var ui_canvas = document.getElementById('ui-canvas');
	ui_canvas.style.top = 0;
	ui_canvas.style.position = 'absolute' 

	this.glRenderer = new THREE.WebGLRenderer( { canvas: ui_canvas,
											     alpha: true,
											     antialias: true } );
	this.glRenderer.setClearColor( 0x000000, 0 ); // the default
	this.glRenderer.setSize(this.width, this.height);
	this.glScene = new THREE.Scene();
	this.glScene.add(this.camera);

	var pointLight = new THREE.PointLight(0xFFFFFF);
	pointLight.position.x = 10;
	pointLight.position.y = 50;
	pointLight.position.z = 130;
	this.glScene.add(pointLight);

	document.body.appendChild(this.glRenderer.domElement);

	this.laser = new LaserPointer(this.glRenderer, this.glScene, this.camera);

	this.raycaster = new THREE.Raycaster();

	this.mouse = new THREE.Vector2();
	this.mouseOffset = new THREE.Vector3();
	this.INTERSECTED = null;
	this.SELECTED = null;

	this.glRenderer.domElement.addEventListener( 'mousemove', this, false );
	this.glRenderer.domElement.addEventListener( 'mousedown', this, false );
	this.glRenderer.domElement.addEventListener( 'mouseup',   this, false );

	this.container = document.getElementById('container');
	{
		this.stats = new Stats();
		this.stats.domElement.style.position = 'absolute';
		this.stats.domElement.style.top = '0px';
		this.container.appendChild( this.stats.domElement );
	}
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
	if (!this.needsReset)
	    return;

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

Renderer.prototype.setup = function(shaderSources)
{
	var gl = GLU.gl;

	// Quad VBO for rendering textures
	this.quadVbo = this.createQuadVbo();

	// shaderSources is a dict from name (e.g. "viewport")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.initProgram  = new GLU.Shader('init',  shaderSources);
	this.traceProgram = new GLU.Shader('trace', shaderSources);
	this.lineProgram  = new GLU.Shader('line',  shaderSources);
	this.compProgram  = new GLU.Shader('comp',  shaderSources);
	this.passProgram  = new GLU.Shader('pass',  shaderSources);
	
	// table of 256 vec4 RGB colors, corresponding to the 256 wavelength samples
	// between 360.0 and 750.0 nanometres
	// @todo:  for now we will just assume a flat emission spectrum, i.e pure white light.
	this.LAMBDA_MIN = 360.0;
    this.LAMBDA_MAX = 750.0;
	this.spectrumTable = wavelengthToRgbTable();
	this.spectrum = new GLU.Texture(this.spectrumTable.length/4, 1, 4, true,  true, true, this.spectrumTable);
  
	//gl.viewport(0, 0, this.width, this.height);
	this.raySize = 128;
	this.resetActiveBlock();
	this.rayCount = this.raySize*this.raySize;
	this.currentState = 0;
	this.needsReset = true;
	this.maxPathLength = 15;
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

	this.fbo = new GLU.RenderTarget();

	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	this.resize(this.width, this.height);
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


Renderer.prototype.composite = function()
{
	this.screenBuffer.bind(0);
	this.compProgram.bind();
	this.quadVbo.bind();
	this.compProgram.uniformTexture("Frame", this.screenBuffer);

	// Tonemap to effectively divide by total number of emitted photons
	// (and also apply gamma correction)
	this.compProgram.uniformF("Exposure", 10.0/(Math.max(this.samplesTraced, 1)));//this.raySize*this.activeBlock)));
	
	this.quadVbo.draw(this.compProgram, this.gl.TRIANGLE_FAN);
}


// With y-axis as the polar north
Renderer.prototype.sphericalPolar = function(theta, phi)
{
	Sp = Math.sin(phi);
	Cp = Math.cos(phi);
	St = Math.sin(theta);
	Ct = Math.cos(theta);
	return new THREE.Vector3( Cp*St, Ct, Sp*St );
}


Renderer.prototype.render = function()
{
	if (!this.initialized) return;

	var gl = this.gl;

	this.needsReset = true;

	var current = this.currentState;
	var next    = 1 - current;

	this.fbo.bind();

	// trace rays
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
		// Start all rays at emission point(s)
		this.initProgram.bind();

		// Read random seed from the current state
		this.rayStates[current].rngTex.bind(0); 
		this.initProgram.uniformTexture("RngData", this.rayStates[current].rngTex);

		// Read wavelength -> RGB table 
		this.spectrum.bind(1);
		this.initProgram.uniformTexture("Spectrum", this.spectrum);
		
		// Emitter data (currently just location and direction)
		emitterPos = new THREE.Vector3(-10.0, 0.0, 0.0);
		this.initProgram.uniform3F("EmitterPos", emitterPos.x, emitterPos.y, emitterPos.z);

		emitterDir = this.sphericalPolar(Math.PI/2.0, 0.0);
		this.initProgram.uniform3F("EmitterDir", emitterDir.x, emitterDir.y, emitterDir.z);

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

		// Use the current state as the initial conditions
		this.rayStates[current].bind(this.traceProgram);

		// Generate the next ray state
		this.quadVbo.draw(this.traceProgram, gl.TRIANGLE_FAN);

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
	if (this.pathLength==this.maxPathLength || this.wavesTraced==0)
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

	// Render laser pointer

	//this.laser.material.needsUpdate = true;
	this.laser.render(this.glScene);

	// Update raytracing state
	this.currentState = next;

	this.stats.update();
}


Renderer.prototype.onDocumentMouseMove = function(event)
{
	event.preventDefault();

	this.mouse.x =   (event.clientX / window.innerWidth)*2 - 1;
	this.mouse.y = - (event.clientY / window.innerHeight)*2 + 1;

	this.raycaster.setFromCamera( this.mouse, this.camera );

	obj = this.laser.getObjects();

	if ( this.SELECTED ) 
	{
		var intersects = this.raycaster.intersectObject( obj['backPlane'] );
		if ( intersects.length > 0 ) 
		{
			this.SELECTED.position.copy( intersects[0].point.sub(this.mouseOffset) );
		}
		return;
	}

	var intersects = this.raycaster.intersectObject(obj['housing']);

	if ( intersects.length > 0 )
	{
		if ( this.INTERSECTED != intersects[0].object )
		{
			this.INTERSECTED = intersects[0].object;
			obj['housing'].position.copy( this.INTERSECTED.position );
			obj['backPlane'].lookAt( this.camera.position );
		}
		this.container.style.cursor = 'pointer';
	}
	else
	{
		this.INTERSECTED = null;
		this.container.style.cursor = 'auto';
	}
}


Renderer.prototype.onDocumentMouseDown = function(event)
{
	event.preventDefault();

	this.raycaster.setFromCamera( this.mouse, this.camera );

	obj = this.laser.getObjects();
	var intersects = this.raycaster.intersectObject(obj['housing']);

	if ( intersects.length > 0 )
	{
		controls.enabled = false;
		
		this.SELECTED = intersects[0].object;

		var intersects = this.raycaster.intersectObject( obj['backPlane'] );
		if ( intersects.length > 0 ) 
		{
			this.mouseOffset.copy( intersects[0].point ).sub( obj['backPlane'].position );
		}

		this.container.style.cursor = 'move';


	}
}


Renderer.prototype.onDocumentMouseUp = function(event)
{
	event.preventDefault();
	controls.enabled = true;

	if ( this.INTERSECTED ) 
	{
		obj = this.laser.getObjects();
		obj['backPlane'].position.copy( this.INTERSECTED.position );
		this.SELECTED = null;
	}

	container.style.cursor = 'auto';
}





