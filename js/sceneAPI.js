
/** 
* @constructor 
*/
function Scene() {}

/////////////////////////////////////////////////////////////////////////////////
// Scene state initialization
/////////////////////////////////////////////////////////////////////////////////

/**
* Optionally (but usually), provide this function to set scene and renderer initial state.
* This is called only once during execution, on loading the scene HTML page (or on global reset via 'R' key).
* @param {Snelly} snelly - The snelly object
*/
Scene.prototype.init = function(snelly)
{
	////////////////// copy-pasted console output on 'O', begin /////////////////////
	
	let renderer  = snelly.getRenderer();
	let camera    = snelly.getCamera();
	let controls  = snelly.getControls();
	let materials = snelly.getMaterials();
		
	this.parameters = {};
	this.parameters.foo = 0.631430584918957;
	this.parameters.foo2 = 0.631430584918957;
	this.parameters.bar = 0.5863284002818887;
	this.animFrame = 0;

	snelly.showGUI(false);
		
	// Camera settings:
	// 		camera is a THREE.PerspectiveCamera object
	// 		controls is a THREE.OrbitControls object
	camera.fov = 40;
	camera.up.set(0, 1, 0);
	camera.position.set(-3.7987995289497105, 3.50565369210148, 3.9561774983519684);
	controls.target.set(0.0371344633810977, -0.030567217009180497, 0.022022000504228156);
	controls.zoomSpeed = 2;
	controls.keyPanSpeed = 100;

	// Renderer settings
	//renderer.width = 1280; // (if either width or height are not specified, render size will be taken from window
	//renderer.height = 720; // and will then auto-resize with the window)
	renderer.renderMode = 'pt';  // The other modes are: 'ao', 'normals'
	renderer.maxBounces = 9;
	renderer.maxMarchSteps = 512;
	renderer.radianceClamp = 0.4355179704016914; // (log scale)
	renderer.skyPower = 4.0;
	renderer.skyTemperature = 6000;
	renderer.exposure = 4.5;
	renderer.gamma = 2.2;
	renderer.whitepoint = 2;
	renderer.goalFPS = 10;

	// Material settings
	let surface = materials.loadSurface();
	surface.roughness = 0.05;
	surface.ior = 1.3530655391120507;
	surface.diffuseAlbedo = [0.5, 0.5, 0.5]; //0.7156862745098039, 0.18985910396453856, 0.0771818531334102];
	surface.specAlbedo = [0.0, 0.0, 0.0]; //0.1470588235294118, 0.1470588235294118, 0.1470588235294118];

	let dielectric = materials.loadDielectric('Diamond');
	dielectric.absorptionColor = [1.0, 1.0, 1.0];
	dielectric.absorptionScale = 1.0; // mfp in multiples of scene scale
	dielectric.roughness = 0.030443974630021145;

	let metal = materials.loadMetal('Gold');
	metal.roughness = 0.05;

	////////////////// copy-pasted console output on 'O', end /////////////////////
}

/**
* Optionally, provide this function which generates the init code to re-generate 
* the current UI parameter settings. This will be dumped to the console (along with 
* the rest of the UI state) on pressing key 'O', allowing the scene and renderer
* state to be tweaked in the UI then saved by copy-pasting code into the init function below.
*/
Scene.prototype.initGenerator = function()
{
	return `
this.parameters = {};
this.parameters.foo = ${this.parameters.foo};
this.parameters.bar = ${this.parameters.bar};
this.frame = 0;
	`; 
}

/**
* Optionally, supply an env-map texture URL (must be a lat-long format image).
* (If this is function not implemented, or it returns the empty string, a uniform
* temperature blackbody sky is used).
* @returns {String}
*/
Scene.prototype.envMap = function()
{
  	return 'https://cdn.rawgit.com/portsmouth/envmaps/74e9d389/HDR_040_Field_Bg.jpg';
}

/**
* Optional name (displayed in UI)
* @returns {String}
*/
Scene.prototype.getName = function() { return "Complete API example"; }

/**
* Optional clickable URL (displayed in UI)
* @returns {String}
*/
Scene.prototype.getURL = function() { return "https://github.com/portsmouth/snelly"; }

/////////////////////////////////////////////////////////////////////////////////////
// Shader code
/////////////////////////////////////////////////////////////////////////////////////

/**
* Returns a chunk of GLSL code defining the SDFs which determine the geometry of uber-surface, metal and dielectric materials in the scene.
* Define also (optionally) functions giving the 3d spatial dependence of the material parameters.
* This function is mandatory!

* The code *must* define at least one of the three functions:
*```glsl
*      float SDF_METAL(vec3 X);
*      float SDF_DIELECTRIC(vec3 X);
*      float SDF_DIFFUSE(vec3 X);
*```
* (If only one or two of these are present, the others will not be rendered).
* Optionally, any of the following functions defining the spatial *modulation* of material parameters can be defined.
* (Any of the user-defined functions below can be omitted, in which case they will be replaced with the default indicated).
*```glsl

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 SURFACE_DIFFUSE_REFLECTANCE(in vec3 X);

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 SURFACE_SPECULAR_REFLECTANCE(in vec3 X);

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float SURFACE_IOR(in vec3 X);

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float SURFACE_ROUGHNESS(in vec3 X);

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 METAL_SPECULAR_REFLECTANCE(in vec3 X);

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float METAL_ROUGHNESS(in vec3 X);

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 DIELECTRIC_SPECULAR_REFLECTANCE(in vec3 X);

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float DIELECTRIC_ROUGHNESS(in vec3 X);
*```
* @returns {String}
*/
Scene.prototype.shader = function()
{
	return `

		uniform float _foo;
		uniform float _foo2;
		uniform float _bar;

		#define TWOPI 6.28318530718

		vec2 rotate(vec2 v, float a) 
		{
		 	return vec2(cos(a)*v.x + sin(a)*v.y, -sin(a)*v.x + cos(a)*v.y); 
		}

		float sdBox(vec3 X, vec3 bounds)                     
		{                                     
		    vec3 d = abs(X) - bounds;
		    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
		} 

		float sdCross( in vec3 p )
		{
			float inf = 1.0e6;
			float da = sdBox(p.xyz, vec3(inf,1.0,1.0));
			float db = sdBox(p.yzx, vec3(1.0,inf,1.0));
			float dc = sdBox(p.zxy, vec3(1.0,1.0,inf));
			return min(da,min(db,dc));
		}

		float menger(vec3 X) 
		{
			float _rotate = _foo2;
			float sd = sdBox(X, vec3(1.0));
			float scale = _foo;
			const int iter = 10;
			for (int i=0; i<iter; i++) 
			{
				X.xy = rotate(X.xy, 2.0*sin(TWOPI*_rotate));
				vec3 a = mod(X*scale, 2.0)-1.0;
				scale *= 3.0;
				vec3 r = abs(1.0 - 3.0*abs(a));
				float da = max(r.x,r.y);
				float db = max(r.y,r.z);
				float dc = max(r.z,r.x);
				float c = (min(da,min(db,dc))-1.0)/scale;
				sd = max(sd, c);
			}
			return sd;
		}	

		float sdSphere(vec3 X, float r)                  
		{                                     
			return length(X) - r;       
		}   

		float sdBox(vec3 X, vec3 bmin, vec3 bmax)                     
		{                            
			vec3 center = 0.5*(bmin + bmax);
			vec3 halfExtents = 0.5*(bmax - bmin);         
			vec3 d = abs(X-center) - halfExtents;
			return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
		} 

		// Union
		float opU( float d1, float d2 ) { return min(d1,d2); }

		// Subtraction
		float opS(float A, float B) { return max(-B, A); }

		// Intersection
		float opI( float d1, float d2 ) { return max(d1,d2); }


		// Any of the user-defined functions below can be omitted, in which case they will be
		// replaced with the default indicated.
		// Note that the reflectances can be treated as [0,1] RGB values in sRGB color space.

		//////////////////////////////////////////////////////
		// Surface
		//////////////////////////////////////////////////////

		float SDF_SURFACE(in vec3 X)                     
		{
			vec3 bmin = vec3(-100.0, -1.0, -100.0);
			vec3 bmax = vec3( 100.0,  0.0,  100.0);
			float floor = sdBox(X, bmin, bmax);
			return opU(floor, menger(X));
		}  

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 SURFACE_DIFFUSE_REFLECTANCE(in vec3 X)
		{
			float albedo = 1.0;
			if (X.y<=0.001)
			{
				float ax = 1.0 - pow(0.5*(1.0 + cos(50.0*X.x)), 100.0);
				float ay = 1.0 - pow(0.5*(1.0 + cos(50.0*X.z)), 100.0);
				albedo = ax*ay;
			}
		    return vec3(albedo);
		}

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 SURFACE_SPECULAR_REFLECTANCE(in vec3 X)
		{
		    return vec3(1.0);
		}

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float SURFACE_IOR(in vec3 X)
		{
		    return 1.0;
		}

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float SURFACE_ROUGHNESS(in vec3 X)
		{
		    return 1.0;
		}


		//////////////////////////////////////////////////////
		// Metal
		//////////////////////////////////////////////////////   
		
		float SDF_METAL(in vec3 X)                     
		{	
			 return 1.0e6;//0.5*menger(X*0.5-vec3(1.0));
		} 

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 METAL_SPECULAR_REFLECTANCE(in vec3 X)
		{
		    return vec3(1.0, 1.0, 1.0);
		}

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float METAL_ROUGHNESS(in vec3 X)
		{
		    return 1.0;
		}


		//////////////////////////////////////////////////////
		// Dielectric
		//////////////////////////////////////////////////////   
		
		float SDF_DIELECTRIC(vec3 X)                     
		{
		    return sdSphere(X, _bar);
		}  

		// space-varying multiplier to the UI-exposed color (defaults to vec3(1.0))
		vec3 DIELECTRIC_SPECULAR_REFLECTANCE(in vec3 X)
		{
		    return vec3(1.0, 1.0, 1.0);
		}

		// space-varying multiplier to the UI-exposed constant (defaults to 1.0)
		float DIELECTRIC_ROUGHNESS(in vec3 X)
		{
		    return 1.0;
		}

		////////////////////////////////////////////////////////////////////// end of shader code
	`;
}

/**
* Optional. Set up gui and callbacks for this scene
* @param {GUI} gui - wrapper for dat.GUI object
*/
Scene.prototype.initGui = function(gui)            
{
  	gui.addParameter(this.parameters, {name: 'foo', min: 0.0, max: 1.0});
  	gui.addParameter(this.parameters, {name: 'foo2', min: 0.0, max: 1.0});
  	gui.addParameter(this.parameters, {name: 'bar', min: 0.0, max: 3.0});
}

/**
* Optional. Called whenever the UI is changed,
/* and must sync the params of the shader with the current UI settings
* @param {Shader} shader - wrapper of webGL fragment shader
*/
Scene.prototype.syncShader = function(shader)
{
	shader.uniformF("_foo", this.parameters.foo);
	shader.uniformF("_foo2", this.parameters.foo2);
	shader.uniformF("_bar", this.parameters.bar);
}

/**
* Optional. Gives the raytracer some indication of the (rough) minimum length scale, 
* so it can set tolerances appropriately. This sets the rough length scale of the smallest 
* resolvable structure. (Note that decreasing this will usually lead to longer render times).
* Defaults to 0.0001.
* @returns {number}
*/
Scene.prototype.getMinScale = function()
{
	return 1.0e-4;
}


/**
* Optional. Gives the raytracer some indication of the (rough) maximum length scale, 
* so it can set tolerances appropriately. The raymarcher will march no further
* from the camera than this scale, thus it acts as the "far plane" distance.
* (Note that increasing this will usually lead to longer render times).
* Defaults to 100.0.
* @returns {number}
*/
Scene.prototype.getMaxScale = function()
{
	return 1.0e2;
}


//////////////////////////////////////////////////////////////////////////////////
// Callbacks
//////////////////////////////////////////////////////////////////////////////////

/** 
 * Optional callback before every frame.
 * Animation rendering logic can be implemented here by updating the scene 
 * programmatically according to the global time since init.
 * @param {Snelly} snelly - The Snelly object
 * @param {WebGLRenderingContext} gl - The webGL context
 */
Scene.prototype.preframeCallback = function(snelly, gl)
{
	let renderer  = snelly.getRenderer();
	let camera    = snelly.getCamera();
	let controls  = snelly.getControls();
	let materials = snelly.getMaterials();
	let gui       = snelly.getGUI();

	let FPS = 24.0;
	let time = this.animFrame/FPS;
	let period = 180.0;
	this.endFrame = period * FPS;

	// animate camera (here a simple 'turntable' orbit about the original cam target)
	let axis = camera.up;

	if (this.animFrame > this.endFrame) this.animFrame = 0;
	if (this.animFrame==0)
	{
		// Do any one-time initial setup of the scene state here
		this.advanceFrame = false;

		this.axisProj = new THREE.Vector3();
		this.axisPerp = new THREE.Vector3();
		let targetToCam = new THREE.Vector3();
		targetToCam.copy(camera.position).sub(controls.target);
		this.axisProj.copy(controls.target);
		this.axisProj.addScaledVector(axis, axis.dot(targetToCam));
		this.axisPerp.copy(targetToCam).sub(this.axisProj);
	}

	///
	// Animate scene state according to current time value here
	// (e.g. update scene, camera, materials or renderer parameters)
	///
	let phase = 2.0*Math.PI*time/period;
	let rot = new THREE.Quaternion();
	rot.setFromAxisAngle(axis, phase);
	let newAxisPerp = new THREE.Vector3();
	newAxisPerp.copy(this.axisPerp);
	newAxisPerp.applyQuaternion(rot);
	
	let newCamPos = new THREE.Vector3();
	newCamPos.copy(controls.target);
	newCamPos.add(this.axisProj);
	newCamPos.add(newAxisPerp);
	camera.position.copy(newCamPos);
	controls.update();

	// animate user scene parameters
	this.parameters.foo  = Math.abs(Math.cos(Math.exp(0.5+0.5*Math.sin(phase))));
	this.parameters.foo2 = Math.abs(Math.cos(Math.exp(0.5-0.5*Math.sin(7.0*phase))));
	this.parameters.bar = 3.0*Math.abs(0.333+Math.cos(Math.exp(0.5+0.5*Math.cos(3.0*phase))));

	// animate materials
	let surface = materials.getSurface();
	surface.roughness = Math.pow(Math.abs(Math.sin(11.0*phase)), 3.0);
	surface.diffuseAlbedo = [phase, 0.1, 1.0-phase];

	let dielectric = materials.getDielectric();
	dielectric.roughness = 0.1*Math.pow(Math.abs(Math.sin(phase)), 5.0);

	// Advance scene state to next anim frame, if we just exported a rendered frame
	if (this.advanceFrame)
	{	
		gui.sync();
		let no_recompile = true;
		renderer.reset(no_recompile);
		this.advanceFrame = false;
	}
}


/** 
 * Optional callback after every frame.
 * Animation rendering logic can be implemented here by updating the scene 
 * programmatically according to the global time since init.
 * @param {Snelly} snelly - The Snelly object
 * @param {WebGLRenderingContext} gl - The webGL context
 */
Scene.prototype.postframeCallback = function(snelly, gl)
{
	// The code here posts the framebuffer pixels to a local server, for sequence rendering.
	let renderer  = snelly.getRenderer();
	let targetSPP = 250.0;

	// User code to post webGL framebuffer data to local server for processing
	if (this.animFrame>=0 && renderer.spp>=targetSPP && this.animFrame<=this.endFrame)
	{
		console.log(this.animFrame);
		let dataURI = gl.canvas.toDataURL();
		var mimetype = dataURI.split(",")[0].split(':')[1].split(';')[0];
		var byteString = atob(dataURI.split(',')[1]);
		var u8a = new Uint8Array(byteString.length);
		for (var i = 0; i < byteString.length; i++) 
		{
			u8a[i] = byteString.charCodeAt(i);
		}
		let blob = new Blob([u8a.buffer], { type: mimetype });
		let r = new XMLHttpRequest();
		r.open('POST', 'http://localhost:3999/' + this.animFrame, false);
		r.send(blob);

		this.advanceFrame = true;
		this.animFrame++;
	}
}
