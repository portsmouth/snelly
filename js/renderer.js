

var RAYBUFFER_W;
var RAYBUFFER_H;

var RAY_BUFFER;
var QUAD_VBO;

var textures;
var framebuffers;
var X_index, Xp_index;


var Renderer = function(glu, shader_programs, width, height) 
{
	this.glu = glu;

	var canvas = document.getElementById('canvas');
	this.gl = canvas.getContext('experimental-webgl', {antialias: true});
	if (!this.gl) 
	{
 		this.gl = canvas.getContext("webgl", {antialias: true});
	}
	this.gl.getExtension('OES_texture_float') 
	var gl = this.gl;
	
	this.rendering_program = shader_programs['viewport'];
	this.raytrace_program = shader_programs['raytrace'];

	this.width = width;
	this.height = height;

	gl.viewport(0, 0, width, height);

	// Initialize THREE.js camera
	{
		var VIEW_ANGLE = 45;
		var ASPECT = width / height;
		var NEAR = 0.1;
		var FAR = 10000;

		this.camera = new THREE.PerspectiveCamera(
							    VIEW_ANGLE,
							    ASPECT,
							    NEAR,
							    FAR);
		this.camera.position.z = 300;
	}

	// Init buffers

	// setup rendering GLSL program
	gl.useProgram(this.rendering_program);
	{
		// create texture coords attributes for our line rendering
		this.RAYBUFFER_W = 512;
		this.RAYBUFFER_H = 512;

		// Create the buffer which defines the ray texture coordinates
		RAY_BUFFER = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, RAY_BUFFER);
		
		var vboData = new Float32Array(this.RAYBUFFER_W*this.RAYBUFFER_H*3);
		for (var i = 0; i < this.RAYBUFFER_W * this.RAYBUFFER_H; ++i)
		{
			var u = ((i % this.RAYBUFFER_H) + 0.5) / this.RAYBUFFER_W;
			var v = (Math.floor(i/this.RAYBUFFER_W) + 0.5) / this.RAYBUFFER_H;
			vboData[i*3 + 0] = u;
			vboData[i*3 + 1] = v;
			vboData[i*3 + 2] = i%2;
		}
		gl.bufferData(gl.ARRAY_BUFFER, vboData, gl.STATIC_DRAW);

		// We'll supply texcoords as floats.
		var texcoordLocation = gl.getAttribLocation(this.rendering_program, "a_texCoord");
		gl.enableVertexAttribArray(texcoordLocation);
		gl.vertexAttribPointer(texcoordLocation, 3, gl.FLOAT, false, 0, 0);
	}

	// Prototype of raytracing logic ..
	{
		// Create raytracing GLSL program
		gl.useProgram(this.raytrace_program);

		// Create a 'X' texture, and an 'Xp' texture
		textures = [];
		framebuffers = [];
		X_index = 0;

		// texture X
		{
			var X = glu.createAndSetupTexture(X_index, this.RAYBUFFER_W, this.RAYBUFFER_H);
			var FBO_X = gl.createFramebuffer();

			textures.push(X);
			framebuffers.push(FBO_X);

			console.log("FBO_X: " + FBO_X);
			
			gl.bindFramebuffer(gl.FRAMEBUFFER, FBO_X);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, X, 0);
		}
	
		// texture Xp
		/*
		Xp_index = 1;
		{
			var Xp = createAndSetupTexture(gl, Xp_texIndex, 2048, 2);
			var FBO_Xp = gl.createFramebuffer();

			textures.push(Xp);
			framebuffers.push(FBO_Xp);
			Xp_index = textures.length;

			gl.bindFramebuffer(gl.FRAMEBUFFER, FBO_Xp);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, textureXp, 0);
		}
		*/

		// Create a 2d-vertex buffer and put a single clipspace rectangle in it (2 triangles)
		// which is what the raytracing program will use to render to texture
		QUAD_VBO = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, QUAD_VBO);
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

		var positionLocation = gl.getAttribLocation(this.raytrace_program, "a_position");
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

	gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffers[X_index]);
	gl.bindTexture(gl.TEXTURE_2D, textures[X_index]);

	gl.useProgram(this.raytrace_program);

	// Tell webgl the viewport setting needed for framebuffer.
	gl.viewport(0, 0, this.RAYBUFFER_W, this.RAYBUFFER_H);

	gl.bindBuffer(gl.ARRAY_BUFFER, QUAD_VBO);

	var positionLocation = gl.getAttribLocation(this.raytrace_program, "a_position");
	gl.enableVertexAttribArray(positionLocation);
	gl.vertexAttribPointer(positionLocation, 3, gl.FLOAT, false, 0, 0);	

	// Will run raytrace program over all fragments
	gl.drawArrays(gl.TRIANGLES, 0, 6);
}


Renderer.prototype.render = function()
{
	var gl = this.gl;

	// Raytrace, to generate X, Xp textures containing the next line segments to render
	this.raytrace();

	// Render those line segments in 3d
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);
	gl.useProgram(this.rendering_program);

	// Setup projection matrix
	var projectionMatrix = this.camera.projectionMatrix.toArray();
	var projectionMatrixLocation = gl.getUniformLocation(this.rendering_program, "u_projectionMatrix");
	gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

	// Setup modelview matrix (to match camera)
	this.camera.updateMatrixWorld();
	var matrixWorldInverse = new THREE.Matrix4();
	matrixWorldInverse.getInverse( this.camera.matrixWorld );
	var modelViewMatrix = matrixWorldInverse.toArray();
	var modelViewMatrixLocation  = gl.getUniformLocation(this.rendering_program, "u_modelViewMatrix");
	gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

	// Draw!
	gl.viewport(0, 0, this.width, this.height);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT);

	// Make sure that u_X sampler in rendering program references the X texture:
	var u_XLoc = gl.getUniformLocation(this.rendering_program, "u_X");
	gl.uniform1i(u_XLoc, X_index);

	// provide texture coordinates for the rectangle.
	var texCoordLocation = gl.getAttribLocation(this.rendering_program, "a_texCoord");
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
	gl.activeTexture(gl.TEXTURE0+X_index);
	gl.bindTexture(gl.TEXTURE_2D, textures[X_index]);

	// For now, draw 100 out of our RAYBUFFER_W*RAYBUFFER_H line segments
	gl.bindBuffer(gl.ARRAY_BUFFER, RAY_BUFFER);

	var texcoordLocation = gl.getAttribLocation(this.rendering_program, "a_texCoord");
	gl.enableVertexAttribArray(texcoordLocation);
	gl.vertexAttribPointer(texcoordLocation, 3, gl.FLOAT, false, 0, 0);

	gl.drawArrays(gl.LINES, 0, 512);
}



