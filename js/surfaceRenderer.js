

var SurfaceRendererState = function(width, height)
{
	var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count 
	var rngData      = new Float32Array(width*height*4); // Random number seed
	var depthData    = new Float32Array(width*height*4); // Packed depth
	for (var i = 0; i<width*height; ++i)
	{
		for (var t = 0; t<4; ++t)
		{
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
	}
	this.radianceTex = new GLU.Texture(width, height, 4, true, false, true, radianceData);
	this.rngTex      = new GLU.Texture(width, height, 4, true, false, true, rngData);
	this.depthTex    = new GLU.Texture(width, height, 4, true, false, true, depthData);
}

SurfaceRendererState.prototype.bind = function(shader)
{
	this.radianceTex.bind(0);
	this.rngTex.bind(1);
	this.depthTex.bind(2);
	shader.uniformTexture("Radiance", this.radianceTex);
	shader.uniformTexture("RngData", this.rngTex);
	shader.uniformTexture("Depth", this.depthTex);
}

SurfaceRendererState.prototype.attach = function(fbo)
{
	var gl = GLU.gl;
	fbo.attachTexture(this.radianceTex, 0);
	fbo.attachTexture(this.rngTex, 1);
	fbo.attachTexture(this.depthTex, 2);
	if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) 
	{
		GLU.fail("Invalid framebuffer");
	}
}

SurfaceRendererState.prototype.detach = function(fbo)
{
	var gl = GLU.gl;
	fbo.detachTexture(0);
	fbo.detachTexture(1);
	fbo.detachTexture(2);
}

SurfaceRendererState.prototype.clear = function(fbo)
{
	// clear radiance buffer
	var gl = GLU.gl;
	fbo.bind();
	fbo.drawBuffers(1);
	fbo.attachTexture(this.radianceTex, 0);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT);
	fbo.unbind();
}


var SurfaceRenderer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;

	this.quadVbo = this.createQuadVbo();
	this.fbo = new GLU.RenderTarget();
	
	// Initialize pathtracing textures and fbo:
	this.currentState = 0;
	this.pathStates = [new SurfaceRendererState(this.width, this.height), 
					   new SurfaceRendererState(this.width, this.height)];

	this.maxMarchSteps = 256;
	this.enable = true;
	this.depthTest = false;
	this.showBounds = false;
	this.surfaceAlpha = 0.5;
	this.renderMode = 'blinn';
	this.specPower = 20.0;
	this.kd1 = [0.25, 0.25, 1.0];
	this.kd2 = [1.0, 0.25, 0.25];

	// Load shaders
	this.shaderSources = GLU.resolveShaderSource(["pathtracer", "tonemapper", "pick"]);
	this.compileShaders();

	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	// Trigger initial buffer generation
	this.resize(this.width, this.height);
}

SurfaceRenderer.prototype.createQuadVbo = function()
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

SurfaceRenderer.prototype.reset = function()
{
	this.compileShaders();
	this.currentState = 0;
	this.pathStates[this.currentState].clear(this.fbo);
}


SurfaceRenderer.prototype.compileShaders = function()
{
	// Inject code for the current scene SDF:
	var sceneObj = snelly.getLoadedScene();
	if (sceneObj == null) return;
	var sdfCode = sceneObj.sdf();

	// Copy the current scene and material routines into the source code
	// of the trace fragment shader
	replacements = {};
	replacements.SDF_FUNC        = sdfCode;
	replacements.MAX_MARCH_STEPS = this.maxMarchSteps;

	switch (this.renderMode)
	{
		case "normals": replacements.LIGHTING_FUNC = `
				vec3 LIGHTING( in vec3 V, in vec3 N )
				{
					return 0.5*(1.0+N);
				}
				`; break;

		case "blinn": replacements.LIGHTING_FUNC = `
				vec3 LIGHTING( in vec3 V, in vec3 N )
				{
					vec3 Lights[4];
					const float oosot = 1.0/sqrt(3.0);
					Lights[0] = vec3( 1.0,  1.0,  1.0)*oosot;
					Lights[1] = vec3(-1.0,  1.0, -1.0)*oosot;
					Lights[2] = vec3(-1.0,  -1.0, -1.0)*oosot;
					Lights[3] = vec3(-1.0,  -1.0, -1.0)*oosot;
					vec3 kd[2];
					kd[0] = vec3(${this.kd1[0]}, ${this.kd1[1]}, ${this.kd1[2]});
					kd[1] = vec3(${this.kd2[0]}, ${this.kd2[1]}, ${this.kd2[2]});
					vec3 ks = vec3(1.0, 1.0, 1.0);
					vec3 ka = vec3(0.005, 0.005, 0.005);
					vec3 C = vec3(0.0, 0.0, 0.0);
					for (int l=0; l<4; ++l)
					{
						vec3 L = Lights[l];
						vec3 H = normalize(L+V);
						float NdotH = dot(N, H);
						float NdotL = dot(N, L);
						C += kd[l]*saturate(NdotL) + ks*pow(saturate(NdotH), float(${this.specPower}));
					}
					return C + ka;
				}
				`; break;									 
	}

	// shaderSources is a dict from name (e.g. "trace")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.pathtraceProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
	this.pickProgram      = new GLU.Shader('pick',       this.shaderSources, replacements);
	this.tonemapProgram   = new GLU.Shader('tonemapper', this.shaderSources, null);
}


SurfaceRenderer.prototype.pick = function(xPick, yPick)
{
	var sceneObj = snelly.getLoadedScene();
	if (sceneObj == null) return null;
	var gl = this.gl;
	gl.viewport(0, 0, 1, 1);

	this.pickProgram.bind();
	sceneObj.syncShader(this.pickProgram); 

	// sync camera info to shader	 
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camY = camera.up.clone();
	camY.transformDirection( camera.matrixWorld );
	var camX = new THREE.Vector3();
	camX.crossVectors(camY, camDir);

	var ndcX = xPick;
	var ndcY = yPick;

	this.pickProgram.uniformF("ndcX", ndcX);
	this.pickProgram.uniformF("ndcY", ndcY);
	this.pickProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	this.pickProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	this.pickProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	this.pickProgram.uniform3Fv("camY", [camY.x, camY.y, camY.z]);
	this.pickProgram.uniformF("camNear", camera.near);
	this.pickProgram.uniformF("camFar", camera.far);
	this.pickProgram.uniformF("camFovy", camera.fov);
	this.pickProgram.uniformF("camZoom", camera.zoom);
	this.pickProgram.uniformF("camAspect", camera.aspect);
	this.pickProgram.uniformF("SceneScale", sceneObj.getScale()); 

	var fbo = new GLU.RenderTarget();
	fbo.bind();
	fbo.drawBuffers(1);

	var pickData = new Float32Array(4);
	this.pickTex = new GLU.Texture(1, 1, 4, true, false, true, pickData);
	fbo.attachTexture(this.pickTex, 0);

	// Trace pick ray
	this.quadVbo.bind();
	this.quadVbo.draw(this.pickProgram, gl.TRIANGLE_FAN);

	// Read floating point hit distance output from pick fragment shader
	var pixels = new Uint8Array(4);
	gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
	pixels = new Float32Array(pixels.buffer);
	fbo.unbind();
	var hitDist = pixels[0];
	if (hitDist < 0.0) return null;

	// Compute hit point given pick NDC, camera, and the hit distance
	var fh = camera.near * Math.tan(0.5*camera.fov*Math.PI/180.0) / camera.zoom;
	var fw = camera.aspect * fh;

	var sX = camX.clone();
	var sY = camY.clone();
	sX.multiplyScalar(-fw*ndcX);
	sY.multiplyScalar(fh*ndcY);
	var s = sX.clone();
	s.add(sY);

	var D = camDir.clone();
	D.multiplyScalar(camera.near);
	D.add(s);
	D.normalize(); // ray direction through camPos, towards hit point
	D.multiplyScalar(hitDist);

	var hitPoint = camPos.clone();
	hitPoint.add(D);
	return hitPoint;
}

SurfaceRenderer.prototype.enabled = function()
{
	return this.enable;
}

SurfaceRenderer.prototype.depthTestEnabled = function()
{
	return this.depthTest;
}

SurfaceRenderer.prototype.render = function()
{
	if (!this.enable) return;

	var sceneObj = snelly.getLoadedScene();
	if (sceneObj == null) return;

	////////////////////////////////////////////////
	/// Pathtracing
	////////////////////////////////////////////////
	var gl = this.gl;
	gl.disable(gl.DEPTH_TEST);
	//gl.viewport(0, 0, this.width, this.height);

	this.pathtraceProgram.bind();

	// sync camera info to shader	 
	var camera = snelly.getCamera();
	var camPos = camera.position.clone();
	var camDir = camera.getWorldDirection();
	var camUp = camera.up.clone();
	camUp.transformDirection( camera.matrixWorld );

	var camX = new THREE.Vector3();
	camX.crossVectors(camUp, camDir);

	this.pathtraceProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
	this.pathtraceProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
	this.pathtraceProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
	this.pathtraceProgram.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
	this.pathtraceProgram.uniformF("camNear", camera.near);
	this.pathtraceProgram.uniformF("camFar", camera.far);
	this.pathtraceProgram.uniformF("camFovy", camera.fov);
	this.pathtraceProgram.uniformF("camZoom", camera.zoom);
	this.pathtraceProgram.uniformF("camAspect", camera.aspect);
	this.pathtraceProgram.uniform2Fv("resolution", [this.width, this.height]);
	this.pathtraceProgram.uniformF("SceneScale", sceneObj.getScale()); 

	this.fbo.bind();
	this.fbo.drawBuffers(3);

	var current = this.currentState;
	var next    = 1 - current;

	gl.disable(gl.BLEND);

	// Read data from the 'current' state
	this.pathStates[current].bind(this.pathtraceProgram);

	// Write data into the 'next' state
	this.pathStates[next].attach(this.fbo);

	// Upload current scene SDF shader parameters
	sceneObj.syncShader(this.pathtraceProgram); 

	// Trace one path per pixel
	this.quadVbo.bind();
	this.quadVbo.draw(this.pathtraceProgram, gl.TRIANGLE_FAN);
	
	//this.pathStates[next].detach(this.fbo);
	this.fbo.unbind();

	////////////////////////////////////////////////
	/// Tonemapping / compositing
	////////////////////////////////////////////////

	this.tonemapProgram.bind();

	var radianceTexCurrent = this.pathStates[next].radianceTex;
	var depthTexCurrent = this.pathStates[next].depthTex;

	radianceTexCurrent.bind(0);	
	depthTexCurrent.bind(1);

	this.tonemapProgram.uniformTexture("Radiance", radianceTexCurrent);
	this.tonemapProgram.uniformTexture("DepthSurface", depthTexCurrent);
	this.tonemapProgram.uniformF("exposure", 1.0);
	this.tonemapProgram.uniformF("invGamma", 1.0);
	this.tonemapProgram.uniformF("alpha", this.surfaceAlpha);

	// Tonemap the 'current' radiance buffer to produce this frame's pixels
	// Get depth texture of light rays, to allow surface to 
	// occlude rays and vice versa
	var lightDepthTex = snelly.getLightTracer().getDepthTexture();
	if (lightDepthTex != null && this.depthTest)
	{
		lightDepthTex.bind(2);
		this.tonemapProgram.uniformTexture("DepthLight", lightDepthTex);
		this.tonemapProgram.uniformI("enableDepthTest", 1);
	}
	else
	{
		this.tonemapProgram.uniformI("enableDepthTest", 0);
	}

	gl.enable(gl.BLEND);
	gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
	gl.blendEquation(gl.FUNC_ADD);

	this.quadVbo.bind();
	this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);

	gl.disable(gl.BLEND);

	// Ping-pong ..
	this.currentState = next;
}


SurfaceRenderer.prototype.resize = function(width, height)
{
	this.width = width;
	this.height = height;

	this.pathStates = [new SurfaceRendererState(this.width, this.height), 
					   new SurfaceRendererState(this.width, this.height)];

	this.reset();
}


