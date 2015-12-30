


//////////////////////////////////////////////////////////////////////
// RayState
//////////////////////////////////////////////////////////////////////

var RayState = function(size) 
{
	this.size = size;

	var posData = new Float32Array(size*size*3); // ray position
	var dirData = new Float32Array(size*size*3); // ray direction
	var rngData = new Float32Array(size*size);   // Random number initial seed
	var rgbData = new Float32Array(size*size*4); // Ray color, and throughput

	for (var i = 0; i < size*size; ++i)
	{
		posData[i*3 + 0] = 0.0;
		posData[i*3 + 1] = 0.0;
		posData[i*3 + 2] = 0.0;

		var theta = Math.random()*Math.PI*2.0;
		dirData[i*3 + 0] = Math.cos(theta);
		dirData[i*3 + 1] = 0.0;
		dirData[i*3 + 2] = Math.sin(theta);

		rngData[i] = i;
		
		for (var t = 0; t < 4; ++t)
			rgbData[i*4 + t] = 0.0;
	}

	this.posTex = new GLU.Texture(size, size, 3, true, false, true, posData);
	this.dirTex = new GLU.Texture(size, size, 3, true, false, true, posData);
	this.rngTex = new GLU.Texture(size, size, 1, true, false, true, rngData);
	this.rgbTex = new GLU.Texture(size, size, 4, true, false, true, rgbData);
}


RayState.prototype.bind = function(shader)
{
	this.posTex.bind(0);
	this.dirTex.bind(1);
	this.rngTex.bind(2);
	this.rgbTex.bind(3);

	shader.uniformTexture("PosData", this.posTex);
	shader.uniformTexture("DirData", this.posTex);
	shader.uniformTexture("RngData", this.rngTex);
	shader.uniformTexture("RgbData", this.rgbTex);
}


RayState.prototype.attach = function(fbo)
{
	fbo.attachTexture(this.posTex, 0);
	fbo.attachTexture(this.dirTex, 1);
	fbo.attachTexture(this.rngTex, 2);
	fbo.attachTexture(this.rgbTex, 3);
}


RayState.prototype.detach = function(fbo)
{
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
	this.width = render_canvas.width;
	this.height = render_canvas.height;
	
	this.initialized = false;

	this.QUAD_VBO = null;
	this.textures = [];
	this.framebuffers = [];
	this.X_index = 0;

	// Initialize THREE.js camera
	{
		var VIEW_ANGLE = 45;
		var ASPECT = this.width / this.height ;
		var NEAR = 0.1;
		var FAR = 10000;
		this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
		this.camera.position.z = 10;

		this.camera.updateProjectionMatrix();
	}

	

	var renderer = this;
	GLU.resolveShaderSource(["init", "trace", "line"],
		function onShadersLoaded(shaderSources)
		{
			if (!renderer.initialized)
			{
				renderer.setup(shaderSources);
				renderer.initialized = true;
			}
		}
	);
}


Renderer.prototype.setup = function(shaderSources)
{
	var gl = GLU.gl;

	// shaderSources is a dict from name (e.g. "viewport")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.initProgram  = new GLU.Shader('init',  shaderSources);
	this.traceProgram = new GLU.Shader('trace', shaderSources);
	this.lineProgram  = new GLU.Shader('line',  shaderSources);
	

	//gl.viewport(0, 0, this.width, this.height);
	this.raySize = 512;
	this.rayCount = this.raySize*this.raySize;
	//this.currentState = 0;
	//this.rayStates = [new RayState(this.raySize), new RayState(this.raySize)];

	// setup rendering GLSL program
	this.lineProgram.bind();
	{
		// Create the buffer which defines the ray texture coordinates
		this.rayVbo = new GLU.VertexBuffer();
		this.rayVbo.addAttribute("TexCoord", 3, gl.FLOAT, false);
		
		this.rayVbo.init(this.rayCount*2);

		var vboData = new Float32Array(this.rayCount*2*3);
		for (var i = 0; i < this.rayCount; ++i)
		{
			var u = ((i % this.raySize) + 0.5) / this.raySize;
			var v = (Math.floor(i/this.raySize) + 0.5) / this.raySize;
			vboData[i*6 + 0] = vboData[i*6 + 3] = u;
			vboData[i*6 + 1] = vboData[i*6 + 4] = v;
			vboData[i*6 + 2] = i%2;
			vboData[i*6 + 5] = 1.0;
		}

		this.rayVbo.copy(vboData);
	}

	// Prototype of raytracing logic ..
	this.traceProgram.bind();
	{
		// Create a 'X' texture, and an 'Xp' texture

		// texture X
		{
			var X = GLU.createAndSetupTexture(this.X_index, this.raySize, this.raySize);
			var FBO_X = gl.createFramebuffer();
			this.textures.push(X);
			this.framebuffers.push(FBO_X);
			
			gl.bindFramebuffer(gl.FRAMEBUFFER, FBO_X);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, X, 0);
		}

		// Create a 2d-vertex buffer and put a single clipspace rectangle in it (2 triangles)
		// which is what the raytracing program will use to render to texture
		this.QUAD_VBO = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, this.QUAD_VBO);
		gl.bufferData(
		    gl.ARRAY_BUFFER,
		    new Float32Array([
		        -1.0, -1.0, 0.0,
		         1.0, -1.0, 0.0,
		        -1.0,  1.0, 0.0,
		        -1.0,  1.0, 0.0,
		         1.0, -1.0, 0.0,
		         1.0,  1.0, 0.0]),
		    gl.STATIC_DRAW);

		var positionLocation = this.traceProgram.getAttribLocation("a_position");
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 3, gl.FLOAT, false, 0, 0);
	}
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
}


Renderer.prototype.raytrace = function() 
{
	var gl = this.gl;

	gl.bindFramebuffer(gl.FRAMEBUFFER, this.framebuffers[this.X_index]);
	gl.bindTexture(gl.TEXTURE_2D, this.textures[this.X_index]);

	this.traceProgram.bind();

	// Tell webgl the viewport setting needed for framebuffer.
	gl.viewport(0, 0, this.raySize, this.raySize);

	gl.bindBuffer(gl.ARRAY_BUFFER, this.QUAD_VBO);

	var positionLocation = this.traceProgram.getAttribLocation("a_position");
	gl.enableVertexAttribArray(positionLocation);
	gl.vertexAttribPointer(positionLocation, 3, gl.FLOAT, false, 0, 0);	

	// Will run raytrace program over all fragments
	gl.drawArrays(gl.TRIANGLES, 0, 6);
}


Renderer.prototype.render = function()
{
	if (!this.initialized) return;

	var gl = this.gl;


/*  dispersion visualization
    -------------------------

	For now, single light =  laser at X0, pointed toward W0,  with assumed emission spectrum

	Make a texture F containing a pixel per photon frequency sample.
	Make vec4 textures X, X' and W, W' also each containing a pixel per frequency sample 
	Initialize textures X = X0, and W = W0 (the location and direction of a 'laser beam' emitter)

	for i = 1 to N/B    (B = batch size)
		
		for each ray sample:   X'_i = trace(X_i, W_i) given F_i
		                       W'_i = sample(X'_i, W_i) given F_i

		drawLine(X_i, X'_i)
		
		X_i = X'_i
		W_i = W'_i
*/


	// Raytrace, to generate X, Xp textures containing the next line segments to render
	this.raytrace();

	// Render those line segments in 3d
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);

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

	// Draw!
	gl.viewport(0, 0, this.width, this.height);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT);

	// Make sure that u_X sampler in rendering program references the X texture:
	var u_XLoc = this.lineProgram.getUniformLocation("u_X");
	gl.uniform1i(u_XLoc, this.X_index);

	// Bind the X texture
	gl.activeTexture(gl.TEXTURE0+this.X_index);
	gl.bindTexture(gl.TEXTURE_2D, this.textures[this.X_index]);

	this.rayVbo.bind();
	this.rayVbo.draw(this.lineProgram, gl.LINES, this.rayCount);
}


