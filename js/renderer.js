

/*
var RAY_BUFFER;
var QUAD_VBO;

var textures;
var framebuffers;
var X_index, Xp_index;
*/

var Renderer = function(width, height) 
{
	this.gl = GLU.gl;
	var gl = GLU.gl;
	
	this.initialized = false;

	this.RAY_BUFFER = null;
	this.QUAD_VBO = null;
	this.textures = [];
	this.framebuffers = [];
	this.X_index = 0;

	// Initialize THREE.js camera
	{
		var VIEW_ANGLE = 45;
		var ASPECT = width / height;
		var NEAR = 0.1;
		var FAR = 10000;
		this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
		this.camera.position.z = 300;
	}

	this.width = width;
	this.height = height;

	var renderer = this;

	GLU.resolveShaderSource(["render", "raytrace"],

		function onShadersLoaded(shaderSources)
		{
			// shaderSources is a dict from name (e.g. "viewport")
			// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
			renderer.renderProgram   = new GLU.Shader('render',   shaderSources);
			renderer.raytraceProgram = new GLU.Shader('raytrace', shaderSources);

			gl.viewport(0, 0, width, height);

			// Init buffers

			// setup rendering GLSL program
			renderer.renderProgram.bind();
			{
				// create texture coords attributes for our line rendering
				renderer.RAYBUFFER_W = 512;
				renderer.RAYBUFFER_H = 512;

				// Create the buffer which defines the ray texture coordinates
				renderer.RAY_BUFFER = gl.createBuffer();
				gl.bindBuffer(gl.ARRAY_BUFFER, renderer.RAY_BUFFER);
				
				var vboData = new Float32Array(renderer.RAYBUFFER_W*renderer.RAYBUFFER_H*3);
				for (var i = 0; i < renderer.RAYBUFFER_W * renderer.RAYBUFFER_H; ++i)
				{
					var u = ((i % renderer.RAYBUFFER_H) + 0.5) / renderer.RAYBUFFER_W;
					var v = (Math.floor(i/renderer.RAYBUFFER_W) + 0.5) / renderer.RAYBUFFER_H;
					vboData[i*3 + 0] = u;
					vboData[i*3 + 1] = v;
					vboData[i*3 + 2] = i%2;
				}
				gl.bufferData(gl.ARRAY_BUFFER, vboData, gl.STATIC_DRAW);

				// We'll supply texcoords as floats.
				var texcoordLocation = renderer.renderProgram.getAttribLocation("a_texCoord");
				//var texcoordLocation = gl.getAttribLocation(this.rendering_program, "a_texCoord");
				gl.enableVertexAttribArray(texcoordLocation);
				gl.vertexAttribPointer(texcoordLocation, 3, gl.FLOAT, false, 0, 0);
			}

			// Prototype of raytracing logic ..
			renderer.raytraceProgram.bind();
			{
				// Create a 'X' texture, and an 'Xp' texture

				// texture X
				{
					var X = GLU.createAndSetupTexture(renderer.X_index, renderer.RAYBUFFER_W, renderer.RAYBUFFER_H);
					var FBO_X = gl.createFramebuffer();
					renderer.textures.push(X);
					renderer.framebuffers.push(FBO_X);
					
					gl.bindFramebuffer(gl.FRAMEBUFFER, FBO_X);
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, X, 0);
				}

				// Create a 2d-vertex buffer and put a single clipspace rectangle in it (2 triangles)
				// which is what the raytracing program will use to render to texture
				renderer.QUAD_VBO = gl.createBuffer();
				gl.bindBuffer(gl.ARRAY_BUFFER, renderer.QUAD_VBO);
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

				var positionLocation = renderer.raytraceProgram.getAttribLocation("a_position");
				//var positionLocation = gl.getAttribLocation(this.raytrace_program, "a_position");
				gl.enableVertexAttribArray(positionLocation);
				gl.vertexAttribPointer(positionLocation, 3, gl.FLOAT, false, 0, 0);
			}

			renderer.initialized = true;
		}
	);

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

	this.raytraceProgram.bind();

	// Tell webgl the viewport setting needed for framebuffer.
	gl.viewport(0, 0, this.RAYBUFFER_W, this.RAYBUFFER_H);

	gl.bindBuffer(gl.ARRAY_BUFFER, this.QUAD_VBO);

	var positionLocation = this.raytraceProgram.getAttribLocation("a_position");
	gl.enableVertexAttribArray(positionLocation);
	gl.vertexAttribPointer(positionLocation, 3, gl.FLOAT, false, 0, 0);	

	// Will run raytrace program over all fragments
	gl.drawArrays(gl.TRIANGLES, 0, 6);
}


Renderer.prototype.render = function()
{
	var gl = this.gl;

	if (!this.initialized) return;

	// Raytrace, to generate X, Xp textures containing the next line segments to render
	this.raytrace();

	// Render those line segments in 3d
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);

	this.renderProgram.bind();

	// Setup projection matrix
	var projectionMatrix = this.camera.projectionMatrix.toArray();
	var projectionMatrixLocation = this.renderProgram.getUniformLocation("u_projectionMatrix");
	//var projectionMatrixLocation = gl.getUniformLocation(this.rendering_program, "u_projectionMatrix");
	gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

	// Setup modelview matrix (to match camera)
	this.camera.updateMatrixWorld();
	var matrixWorldInverse = new THREE.Matrix4();
	matrixWorldInverse.getInverse( this.camera.matrixWorld );
	var modelViewMatrix = matrixWorldInverse.toArray();
	var modelViewMatrixLocation = this.renderProgram.getUniformLocation("u_modelViewMatrix");
	//var modelViewMatrixLocation  = gl.getUniformLocation(this.rendering_program, "u_modelViewMatrix");
	gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

	// Draw!
	gl.viewport(0, 0, this.width, this.height);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT);

	// Make sure that u_X sampler in rendering program references the X texture:
	var u_XLoc = this.renderProgram.getUniformLocation("u_X");
	//var u_XLoc = gl.getUniformLocation(this.rendering_program, "u_X");
	gl.uniform1i(u_XLoc, this.X_index);

	// provide texture coordinates for the rectangle.
	//var texCoordLocation = gl.getAttribLocation(this.rendering_program, "a_texCoord");
	var texCoordLocation = this.renderProgram.getAttribLocation("a_texCoord");
	var texCoordBuffer = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([
		  0.0,  0.0,
		  1.0,  0.0,
		  0.0,  1.0,
		  0.0,  1.0,
		  1.0,  0.0,
		  1.0,  1.0]), gl.STATIC_DRAW);
	gl.enableVertexAttribArray(texCoordLocation);
	gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);
	 
	// Bind the X texture
	gl.activeTexture(gl.TEXTURE0+this.X_index);
	gl.bindTexture(gl.TEXTURE_2D, this.textures[this.X_index]);

	// For now, draw 100 out of our RAYBUFFER_W*RAYBUFFER_H line segments
	gl.bindBuffer(gl.ARRAY_BUFFER, this.RAY_BUFFER);

	var texCoordLocation = this.renderProgram.getAttribLocation("a_texCoord");
	gl.enableVertexAttribArray(texCoordLocation);
	gl.vertexAttribPointer(texCoordLocation, 3, gl.FLOAT, false, 0, 0);

	gl.drawArrays(gl.LINES, 0, 512);
}



