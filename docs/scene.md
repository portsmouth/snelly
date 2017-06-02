
## Scene

### Geometry

We define the rendered scene geometry by specifying, via the {@link Scene#shader}, three GLSL functions:

	- `SDF_SURFACE(vec3 X)`: the SDF of the uber-surface material
	- `SDF_METAL(vec3 X)`: the SDF of the (selected) metal material
	- `SDF_DIELECTRIC(vec3 X)`: the SDF of the (selected) dielectric material
  
These functions are assumed to be SDFs where the negative region corresponds to the interior of the body.
Thus there are at most only three types of material in the scene.

The details of the properties of the three material types can be selected via the {@link Materials} object.
In addition, spatial dependence of the material surface properties can be introduced by providing modulating GLSL functions.


### Lighting

For simplicity, the only light in the scene is a (non-HDRI) environment map. This can be specified via a URL to 
to a lat-long map, via the {@link Scene#envMap} call:
```javascript
	/**
	* Optionally, supply an env-map texture URL (must be a lat-long format image).
	* (If this is function not implemented, or it returns the empty string, a uniform
	* temperature blackbody sky is used).
	*/
	Scene.prototype.envMap = function()
	{
	  	return 'https://cdn.rawgit.com/portsmouth/envmaps/74e9d389/HDR_040_Field_Bg.jpg';
	  	//return 'https://cdn.rawgit.com/portsmouth/envmaps/7405220b/HDR_041_Path_Bg.jpg';
	  	//return 'https://cdn.rawgit.com/portsmouth/envmaps/74e9d389/HDR_112_River_Road_2_Bg.jpg';
	  	//return 'https://cdn.rawgit.com/portsmouth/envmaps/7405220b/HDR_110_Tunnel_Bg.jpg';
	  	//return 'https://cdn.rawgit.com/portsmouth/envmaps/7405220b/HDR_Free_City_Night_Lights_Bg.jpg';
	}
```

Or otherwise will be taken to be a constant intensity sky. 
In both cases, the sky spectrum is modulated by a blackbody emission spectrum with adjustable temperature (Set via {@link Renderer#skyTemperature}).


### Saving scene state

Often we want to explore and fine-tune a scene by moving the camera around, tweaking the scene parameters and materials, and adjusting renderer settings. This work would be wasted without a way to save the resulting scene state. This is provided by the simple mechanism of pressing the 'O' key to dump to the console a Javascript code which can be inserted into the {@link Scene#init} function, to replicate the scene state. An example of this output is:
```javascript
/******* copy-pasted console output on 'O', begin *******/

let renderer  = snelly.getRenderer();
let camera    = snelly.getCamera();
let controls  = snelly.getControls();
let materials = snelly.getMaterials();
	
snelly.showGUI(true);

/** Camera settings:
*		camera is a THREE.PerspectiveCamera object
* 		controls is a THREE.OrbitControls object
*/
camera.fov = 45;
camera.up.set(0, 1, 0);
camera.position.set(10.32444633217751, -0.5153825294569127, -6.7722559545030085);
controls.target.set(-0.5561757324087446, -0.6317558415254627, -1.159897890031697);
controls.zoomSpeed = 2;
controls.keyPanSpeed = 100;

/** Renderer settings **/
renderer.renderMode = 'pt';  // The other modes are: 'ao', 'normals'
renderer.maxBounces = 5;
renderer.maxMarchSteps = 256;
renderer.radianceClamp = 3; // (log scale)
renderer.skyPower = 1;
renderer.skyTemperature = 6000;
renderer.exposure = 3;
renderer.gamma = 2.2;
renderer.whitepoint = 2;
renderer.goalFPS = 20;

/** Material settings **/
let surface = materials.loadSurface();
surface.roughness = 0.1;
surface.ior = 1.5;
surface.diffuseAlbedo = [1, 1, 1];
surface.specAlbedo = [1, 1, 1];

let dielectric = materials.loadDielectric('Diamond');
dielectric.absorptionColor = [0.5, 0.5, 0.5];
dielectric.absorptionScale = 100; // mfp in multiples of scene scale
dielectric.roughness = 0.007431448254773935;

let metal = materials.loadMetal('Copper');
metal.roughness = 0.006691306063588591;

/******* copy-pasted console output on 'O', end *******/
```

In order for user-specified scene parameters to be saved this way, it is necessay to implement a function to output the appropriate regeneration code. For example if the scene uses parameters
```javascript
	this.parameters = {};
	this.parameters.foo = 0.63143;
	this.parameters.foo2 = 0.631;
	this.parameters.bar = 0.586;
	this.animFrame = 0;
```

the appropriate generation code would be:
```javascript
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
```
With this code in place, the output on pressing 'O' is then a faithful representation of the entire scene state.



### Callbacks and animation

For implementation of custom animation logic, we use the simple mechanism of pre- and post-frame user callbacks, wherein the user can implement whatever logic he needs. See the provided examples for details of how to use this implement animating scenes, and movie rendering.

```
Scene.prototype.preframeCallback = function(snelly, gl);
Scene.prototype.postframeCallback = function(snelly, gl);
```

