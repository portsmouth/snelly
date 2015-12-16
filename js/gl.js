




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
function createProgramFromScripts(gl, vertexShaderFile, fragmentShaderFile)
{
	var vertexShader;
	var fragmentShader;
	var crap;

	//console.log("vertexShaderFile: " + vertexShaderFile);
	//console.log("fragmentShaderFile: " + fragmentShaderFile);

	require(
		["text!shaders/"+vertexShaderFile, "text!shaders/"+fragmentShaderFile],
		function(vertexShaderCode, fragmentShaderCode) {

			//console.log("vertexShaderCode: " + vertexShaderCode);
			//console.log("fragmentShaderCode: " + fragmentShaderCode);

			vertexShader   = compileShader(gl, vertexShaderFile, vertexShaderCode, gl.VERTEX_SHADER);
			fragmentShader = compileShader(gl, fragmentShaderFile, fragmentShaderCode, gl.FRAGMENT_SHADER);

			//console.log("vertexShader: " + vertexShader);
			//console.log("fragmentShader: " + fragmentShader);

			window[vertexShaderFile]   = vertexShader;
			window[fragmentShaderFile] = fragmentShader;

			//console.log(window);
		}
	);

	checkLoaded(vertexShaderFile, fragmentShaderFile);

	//
	// OK ----   clearly need to actually be requesting to load the entire set of shaders asynchonously.
	// *Then* in the callback indicating that has finished, continue to the next phase.
	//


	console.log(window[vertexShaderFile]);

	console.log("Loaded " + vertexShaderFile + " and " + fragmentShaderFile);

	return createProgram(gl, vertexShader, fragmentShader);
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



