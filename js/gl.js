

var GLU = {};

(function() {

    var canvas = document.getElementById('canvas');
	this.gl = canvas.getContext('experimental-webgl', {antialias: true});
	if (!this.gl) 
	{
 		this.gl = canvas.getContext("webgl", {antialias: true});
	}
 
	this.gl.getExtension('OES_texture_float')

 	///////////////////////////////////////////////////
 	// GLU namespace functions
 	///////////////////////////////////////////////////

 	this.createAndSetupTexture = function(textureUnitIndex, width, height) 
	{
		var gl = this.gl;
		var texture = gl.createTexture();

		gl.activeTexture(gl.TEXTURE0+textureUnitIndex);
		gl.bindTexture(gl.TEXTURE_2D, texture);

		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.FLOAT, null);

		return texture;
	}

	this.resolveShaderSource = function(shader_names, onFinishedCallback)
	{
		var shader_files = [];
		for (var i = 0; i < shader_names.length; i++) 
		{
			shader_files.push( "text!shaders/"+shader_names[i]+'-vertex-shader.glsl' );
			shader_files.push( "text!shaders/"+shader_names[i]+'-fragment-shader.glsl' );
		}

		require(
			shader_files,
			function() 
			{
				var shaderSources = {};
				var n = 0;
				for (var i = 0; i<(arguments.length)/2; i++) 
				{
					shaderSources[shader_names[i]] = 
					{
						'v': arguments[n],
						'f': arguments[n+1]
					};
					n += 2;
				}
				onFinishedCallback(shaderSources);
			}
		);
	}

	this.createProgram = function(vertexShader, fragmentShader) 
	{
		var gl = this.gl;

		// create a program.
		var program = gl.createProgram();

		// attach the shaders.
		gl.attachShader(program, vertexShader);
		gl.attachShader(program, fragmentShader);

		// link the program.
		gl.linkProgram(program);

		// Check if it linked.
		var success = gl.getProgramParameter(program, gl.LINK_STATUS);
		if (!success) 
		{
			// something went wrong with the link
			throw ("Program failed to link, error: " + gl.getProgramInfoLog (program));
		}

		return program;
	}

	this.compileShaderSource = function(shaderName, shaderSource, shaderType)
	{
		var gl = this.gl;
		var shader = gl.createShader(shaderType); // Create the shader object
		gl.shaderSource(shader, shaderSource); // Set the shader source code.
		gl.compileShader(shader);              // Compile the shader
		var success = gl.getShaderParameter(shader, gl.COMPILE_STATUS); // Check if it compiled
		if (!success) 
		{
			// Something went wrong during compilation; get the error
			throw ("Could not compile " + shaderName + " shader:" + gl.getShaderInfoLog(shader));
		}
		return shader;
	}


 	///////////////////////////////////////////////////
 	// GLU.Shader object
 	///////////////////////////////////////////////////

	this.Shader = function(name, shaderSources)
	{
		var gl = GLU.gl;
		shaderSource = shaderSources[name];
		vertSource = shaderSource.v;
		fragSource = shaderSource.f;
		var vertexShader       = GLU.compileShaderSource(name, vertSource, gl.VERTEX_SHADER);
		var fragmentShader     = GLU.compileShaderSource(name, fragSource, gl.FRAGMENT_SHADER);
		this.program = GLU.createProgram(vertexShader, fragmentShader);
	}

	this.Shader.prototype.bind = function() 
	{
		GLU.gl.useProgram(this.program);
	}

	this.Shader.prototype.getAttribLocation = function(attribName)
	{
		return GLU.gl.getAttribLocation(this.program, attribName);
	}

	this.Shader.prototype.getUniformLocation = function(uniformName)
	{
		return GLU.gl.getUniformLocation(this.program, uniformName);
	}


}).apply(GLU);  








