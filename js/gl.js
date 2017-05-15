

var GLU = {};

(function() {

	///////////////////////////////////////////////////
	// GLU namespace functions
	///////////////////////////////////////////////////

	this.setupGL = function()
	{
		try 
		{
			var gl = this.canvas.getContext("webgl") || this.canvas.getContext("experimental-webgl");
		} catch (e) {}
		if (!gl) throw new Error("Could not initialise WebGL");
		this.gl = gl;

		//console.log('Supported webGL extensions: ' + gl.getSupportedExtensions());
		this.floatExt       = gl.getExtension("OES_texture_float");
		this.floatLinExt    = gl.getExtension("OES_texture_float_linear");
		this.floatBufExt    = gl.getExtension("WEBGL_color_buffer_float");
		this.multiBufExt    = gl.getExtension("WEBGL_draw_buffers");
		this.depthTexExt    = gl.getExtension("WEBGL_depth_texture");
		this.blendMinMaxExt = gl.getExtension("EXT_blend_minmax");
		this.sRGBExt        = gl.getExtension('EXT_sRGB');

		if (!this.floatExt || !this.floatLinExt) throw new Error("Your platform does not support float textures");
		if (!this.multiBufExt)                   throw new Error("Your platform does not support the draw buffers extension");
		if (gl.getParameter(this.multiBufExt.MAX_DRAW_BUFFERS_WEBGL) < 4)
		{
			throw new Error("Your platform does not support 4 draw buffers");
		}
	}

	this.glTypeSize = function(type) 
	{
		switch (type) 
		{
			case gl.BYTE:
			case gl.UNSIGNED_BYTE:
			    return 1;
			case gl.SHORT:
			case gl.UNSIGNED_SHORT:
			    return 2;
			case gl.INT:
			case gl.UNSIGNED_INT:
			case gl.FLOAT:
			    return 4;
			default:
			    return 0;
		}
	}

	this.resolveShaderSource = function(shader_names)
	{
		var shaderSources = {};
		for (var i=0; i<shader_names.length; i++)
		{
			var name = shader_names[i];
			shaderSources[name] =
			{
				'v': Shaders[name+'-vertex-shader'],
				'f': Shaders[name+'-fragment-shader']
			};
		}
		return shaderSources;
	}

	this.createProgram = function(vertexShader, fragmentShader) 
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
		if (!success) 
		{
			// something went wrong with the link
			throw ("Program failed to link, error: " + gl.getProgramInfoLog (program));
		}

		return program;
	}

	this.compileShaderSource = function(shaderName, shaderSource, shaderType)
	{
		var shader = gl.createShader(shaderType); // Create the shader object
		gl.shaderSource(shader, shaderSource); // Set the shader source code.
		gl.compileShader(shader);              // Compile the shader
		var success = gl.getShaderParameter(shader, gl.COMPILE_STATUS); // Check if it compiled
		if (!success) 
		{
			// Something went wrong during compilation; get the error
			var shaderTypeStr = (shaderType==gl.VERTEX_SHADER) ? 'vertex' : 'fragment';
			console.log(shaderSource);
			throw ("Could not compile " + shaderName + " " + shaderTypeStr + " shader: " + gl.getShaderInfoLog(shader));
		}
		return shader;
	}


	///////////////////////////////////////////////////
	// GLU.Shader object
	///////////////////////////////////////////////////

	this.Shader = function(name, shaderSources, replacements)
	{
		shaderSource = shaderSources[name];

		// replacements is an optional object whose key-val pairs
		// define patterns to replace in the vertex and fragment shaders
		vertSource = (' ' + shaderSource.v).slice(1);
		fragSource = (' ' + shaderSource.f).slice(1);
		if (replacements != null)
		{
			for (var pattern in replacements) 
			{
				if (replacements.hasOwnProperty(pattern)) 
				{
					vertSource = vertSource.replace(new RegExp(pattern, 'g'), replacements[pattern]);
					fragSource = fragSource.replace(new RegExp(pattern, 'g'), replacements[pattern]);
				}
			}
		};
		var vertexShader       = GLU.compileShaderSource(name, vertSource, gl.VERTEX_SHADER);
		var fragmentShader     = GLU.compileShaderSource(name, fragSource, gl.FRAGMENT_SHADER);
		this.program = GLU.createProgram(vertexShader, fragmentShader);
		if (!gl.getProgramParameter(this.program, gl.LINK_STATUS))
			alert("Could not initialise shaders");
		this.uniforms = {};
	}

	this.Shader.prototype.bind = function() 
	{
		gl.useProgram(this.program);
	}

	this.Shader.prototype.getAttribLocation = function(attribName)
	{
		return gl.getAttribLocation(this.program, attribName);
	}

	this.Shader.prototype.getUniformLocation = function(uniformName)
	{
		return gl.getUniformLocation(this.program, uniformName);
	}

	this.Shader.prototype.uniformIndex = function(name) 
	{
	    if (!(name in this.uniforms))
	        this.uniforms[name] = gl.getUniformLocation(this.program, name);
	    return this.uniforms[name];
	}

	this.Shader.prototype.uniformTexture = function(name, texture) 
	{
	    var id = this.uniformIndex(name);
	    if (id != -1)
	        gl.uniform1i(id, texture.boundUnit);
	}

	this.Shader.prototype.uniformI = function(name, i) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform1i(id, i);
	}

	this.Shader.prototype.uniformF = function(name, f) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform1f(id, f);
	}

	this.Shader.prototype.uniform2F = function(name, f1, f2) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform2f(id, f1, f2);
	}

	this.Shader.prototype.uniform1Fv = function(name, fvec) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform1fv(id, fvec);
	}

	this.Shader.prototype.uniform2Fv = function(name, fvec2) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform2fv(id, fvec2);
	}

	this.Shader.prototype.uniform3F = function(name, f1, f2, f3) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform3f(id, f1, f2, f3);
	}


	this.Shader.prototype.uniform3Fv = function(name, fvec3) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform3fv(id, fvec3);
	}

	this.Shader.prototype.uniform4F = function(name, f1, f2, f3, f4) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform4F(id, f1, f2, f3, f4);
	}

	this.Shader.prototype.uniform4Fv = function(name, fvec4) 
	{
		var id = this.uniformIndex(name);
		if (id != -1)
		    gl.uniform4fv(id, fvec4);
	}

	///////////////////////////////////////////////////
	// GLU.VertexBuffer object
	///////////////////////////////////////////////////

	this.VertexBuffer = function()
	{
		this.attributes = [];
		this.elementSize = 0;
	}

	this.VertexBuffer.prototype.bind = function()
	{
		gl.bindBuffer(gl.ARRAY_BUFFER, this.glName);
	}

	this.VertexBuffer.prototype.addAttribute = function(name, size, type, norm)
	{
		this.attributes.push({
		    "name": name,
		    "size": size,
		    "type": type,
		    "norm": norm,
		    "offset": this.elementSize,
		    "index": -1
		});
		this.elementSize += size*GLU.glTypeSize(type);
	}

	this.VertexBuffer.prototype.init = function(numVerts) 
	{
		this.length = numVerts;
		this.glName = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, this.glName);
		gl.bufferData(gl.ARRAY_BUFFER, this.length*this.elementSize, gl.STATIC_DRAW);
	}

	this.VertexBuffer.prototype.copy = function(data) 
	{
		if (data.byteLength != this.length*this.elementSize)
		{
			throw new Error("Resizing VBO during copy strongly discouraged");
		}
		gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
	}

	this.VertexBuffer.prototype.draw = function(shader, mode, length) 
	{
		for (var i = 0; i < this.attributes.length; ++i) 
		{
			this.attributes[i].index = gl.getAttribLocation(shader.program, this.attributes[i].name);
			if (this.attributes[i].index >= 0)
			{
				var attr = this.attributes[i];
				gl.enableVertexAttribArray(attr.index);
				gl.vertexAttribPointer(attr.index, attr.size, attr.type, attr.norm, this.elementSize, attr.offset);
			}
		}

		gl.drawArrays(mode, 0, length ? length : this.length);

		for (var i = 0; i < this.attributes.length; ++i)
		{
			if (this.attributes[i].index >= 0) 
			{
				gl.disableVertexAttribArray(this.attributes[i].index);
				this.attributes[i].index = -1;
			}
		}
	}


	///////////////////////////////////////////////////
	// GLU.Texture object
	///////////////////////////////////////////////////

	this.Texture = function(width, height, channels, isFloat, isLinear, isClamped, texels) 
	{
		var coordMode = isClamped ? gl.CLAMP_TO_EDGE : gl.REPEAT;
		this.type     = isFloat   ? gl.FLOAT         : gl.UNSIGNED_BYTE;
		this.format   = [gl.LUMINANCE, gl.RG, gl.RGB, gl.RGBA][channels - 1];

		this.width  = width;
		this.height = height;

		this.glName = gl.createTexture();
		gl.bindTexture(gl.TEXTURE_2D, this.glName);
		gl.texImage2D(gl.TEXTURE_2D, 0, this.format, this.width, this.height, 0, this.format, this.type, texels);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, coordMode);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, coordMode);
		this.setSmooth(isLinear);

		this.boundUnit = -1;
	}

	this.Texture.prototype.setSmooth = function(smooth) 
	{
		var interpMode = smooth ? gl.LINEAR : gl.NEAREST;
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, interpMode);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, interpMode);
	}

	this.Texture.prototype.copy = function(texels) 
	{
		gl.texImage2D(gl.TEXTURE_2D, 0, this.format, this.width, this.height, 0, this.format, this.type, texels);
	}

	this.Texture.prototype.bind = function(unit) 
	{
		gl.activeTexture(gl.TEXTURE0 + unit);
		gl.bindTexture(gl.TEXTURE_2D, this.glName);
		this.boundUnit = unit;
	}


	// creates a texture info { width: w, height: h, texture: tex }
	// The texture will start with 1x1 pixels and be updated
	// when the image has loaded
	this.loadImageAndCreateTextureInfo = function(url, callback) 
	{
		var tex = gl.createTexture();
		gl.bindTexture(gl.TEXTURE_2D, tex);
		if (GLU.sRGBExt != null) gl.texImage2D(gl.TEXTURE_2D, 0, GLU.sRGBExt.SRGB_EXT, 1, 1, 0, GLU.sRGBExt.SRGB_EXT, gl.UNSIGNED_BYTE, new Uint8Array([0, 0, 255, 255]));
		else                     gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([0, 0, 255, 255]));
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
		var imgInfo = {
			width: 1,   // we don't know the size until it loads
			height: 1,
			texture: tex,
			url: url,
		};
		var img = new Image();
		img.addEventListener('load', function() {
			imgInfo.width = img.width;
			imgInfo.height = img.height;
			imgInfo.url = url;
			imgInfo.tex = tex;
			gl.bindTexture(gl.TEXTURE_2D, tex);
			if (GLU.sRGBExt != null) gl.texImage2D(gl.TEXTURE_2D, 0, GLU.sRGBExt.SRGB_EXT, GLU.sRGBExt.SRGB_EXT, gl.UNSIGNED_BYTE, img);
			else                     gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, img);
			callback(imgInfo);
		});
		if ((new URL(url)).origin !== window.location.origin) 
		{
  		  	img.crossOrigin = "";
  		}
		img.src = url;
		return imgInfo;
	}

	///////////////////////////////////////////////////
	// GLU.RenderTarget object
	///////////////////////////////////////////////////

	this.RenderTarget = function() 
	{
		this.glName = gl.createFramebuffer();
	}

	this.RenderTarget.prototype.bind = function() 
	{
		gl.bindFramebuffer(gl.FRAMEBUFFER, this.glName);
	}

	this.RenderTarget.prototype.unbind = function() 
	{
		gl.bindFramebuffer(gl.FRAMEBUFFER, null);
	}

	this.thing = false

	this.RenderTarget.prototype.attachTexture = function(texture, index) 
	{
		gl.framebufferTexture2D(gl.FRAMEBUFFER, GLU.multiBufExt.COLOR_ATTACHMENT0_WEBGL + index, gl.TEXTURE_2D, texture.glName, 0);
	}

	this.RenderTarget.prototype.detachTexture = function(index) 
	{
		gl.framebufferTexture2D(gl.FRAMEBUFFER, GLU.multiBufExt.COLOR_ATTACHMENT0_WEBGL + index, gl.TEXTURE_2D, null, 0);
	}

	this.RenderTarget.prototype.drawBuffers = function(numBufs) 
	{
		var buffers = [];
		for (var i = 0; i<numBufs; ++i)
		    buffers.push(GLU.multiBufExt.COLOR_ATTACHMENT0_WEBGL + i);
		GLU.multiBufExt.drawBuffersWEBGL(buffers);
	}


	///////////////////////////////////////////////////
	// GLU logic
	///////////////////////////////////////////////////

	this.fail = function(message)
	{
		var sorryP = document.createElement("p"); 
		sorryP.appendChild(document.createTextNode("Sorry! :("));
		sorryP.style.fontSize = "50px";

		var failureP = document.createElement("p");
		failureP.className = "warning-box";
		failureP.innerHTML = message;

		var failureDiv = document.createElement("div"); 
		failureDiv.className = "center";
		failureDiv.appendChild(sorryP);
		failureDiv.appendChild(errorImg);
		failureDiv.appendChild(failureP);

		document.getElementById("content").appendChild(failureDiv);
		this.overlay.style.display = this.canvas.style.display = 'none';
	}

	this.canvas = document.getElementById('render-canvas');
	this.canvas.width = 1;
	this.canvas.height= 1;

	try 
	{
		this.setupGL();
	}
	catch (e) 
	{
		/* GL errors at this stage are to be expected to some degree,
		   so display a nice error message and call it quits */
		this.fail(e.message + ". This can't run in your browser.");
		return;
	}

	var gl = this.gl;

}).apply(GLU);  








