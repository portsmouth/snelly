




/**
 * Creates and compiles a shader.
 *
 * @param {!WebGLRenderingContext} gl The WebGL Context.
 * @param {string} shaderSource The GLSL source code for the shader.
 * @param {number} shaderType The type of shader, VERTEX_SHADER or
 *     FRAGMENT_SHADER.
 * @return {!WebGLShader} The shader.
 */
function compileShader(gl, shaderName, shaderSource, shaderType) 
{
	// Create the shader object
	var shader = gl.createShader(shaderType);

	gl.shaderSource(shader, shaderSource); // Set the shader source code.
	gl.compileShader(shader);              // Compile the shader

	// Check if it compiled
	var success = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
	if (!success) 
	{
		// Something went wrong during compilation; get the error
		throw "could not compile " + shaderName + " shader:" + gl.getShaderInfoLog(shader);
	}

	return shader;
}




/**
 * Creates a program from 2 shaders.
 *
 * @param {!WebGLRenderingContext) gl The WebGL context.
 * @param {!WebGLShader} vertexShader A vertex shader.
 * @param {!WebGLShader} fragmentShader A fragment shader.
 * @return {!WebGLProgram} A program.
 */
function createProgram(gl, vertexShader, fragmentShader) 
{
	// create a program.
	var program = gl.createProgram();

	// attach the shaders.
	gl.attachShader(program, vertexShader);
	gl.attachShader(program, fragmentShader);

	// link the program.
	gl.linkProgram(program);

	// Check if it linked.
	var success = gl.getProgramParameter(program, gl.LINK_STATUS);
	if (!success) {
	  // something went wrong with the link
	  throw ("program filed to link:" + gl.getProgramInfoLog (program));
	}

	return program;
};


function checkLoaded(vertexShaderFile, fragmentShaderFile)
{
	if( !(vertexShaderFile in window) || !(fragmentShaderFile in window) )
	{
		setTimeout(checkLoaded, 100);
		return;
	}
}

/**
 * Creates a program from 2 script tags.
 *
 * @param {!WebGLRenderingContext} gl The WebGL Context.
 * @param {string} vertexShaderId The id of the vertex shader script tag.
 * @param {string} fragmentShaderId The id of the fragment shader script tag.
 * @return {!WebGLProgram} A program
 */
function createProgramsFromScripts(gl, shader_names)
{
	var shader_files = [];
	for (var i = 0; i < shader_names.length; i++) 
	{
		var name = shader_names[i];
		shader_files.push( "text!shaders/"+name+'-vertex-shader.glsl' );
		shader_files.push( "text!shaders/"+name+'-fragment-shader.glsl' );
	}

	//console.log("shader_files: " +  shader_files)

	require(

		shader_files,

		function() 
		{
			console.log("hello: " + gl)
			var shader_programs = {};

			//console.log("arguments: " +  arguments);

			var n = 0;
			for (var i = 0; i<(arguments.length)/2; i++) 
			{
				vertexShaderCode   = arguments[n];
				fragmentShaderCode = arguments[n+1]

				//console.log("vertexShaderCode: " + vertexShaderCode);
				//console.log("fragmentShaderCode: " + arguments);

				vertexShader       = compileShader(gl, shader_names[i], vertexShaderCode, gl.VERTEX_SHADER);
				fragmentShader     = compileShader(gl, shader_names[i], fragmentShaderCode, gl.FRAGMENT_SHADER);

				shader_programs[ shader_names[i] ] = createProgram(gl, vertexShader, fragmentShader);
				n += 2;
			}

			init(shader_programs);
		}
	);
}

function createAndSetupTexture(gl, textureUnitIndex, width, height) 
{
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



