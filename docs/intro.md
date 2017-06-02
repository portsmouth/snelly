
### HTML structure

A Snelly scene is a single, standalone HTML file, which has the following overall structure:
```html
<body onload="onLoad();">
<script type="text/javascript" src="../js/compiled/snelly.min.js"></script>
<script type="text/javascript">

function Scene() {}
Scene.prototype.shader = function() 
{ 
  /* GLSL code */ 
  return `
    uniform float foo; 
    float SDF_SURFACE(vec3 X)    { /* <code omitted> */ }
    float SDF_METAL(vec3 X)      { /* <code omitted> */ }
    float SDF_DIELECTRIC(vec3 X) { /* <code omitted> */ }
    vec3 SURFACE_DIFFUSE_REFLECTANCE(in vec3 X) { /* <code omitted> */ }
    // etc.
  `; 
}

Scene.prototype.init = function(snelly) 
{ 
  /* initial scene setup */
  let renderer  = snelly.getRenderer();
  let camera    = snelly.getCamera();
  let controls  = snelly.getControls();
  let materials = snelly.getMaterials();
  
  camera.position.set(6.0, 3.0, -6.0);
  controls.target.set(0.0, 0.0, 0.0);
  renderer.maxBounces = 9;
  materials.loadMetal('Copper');
  materials.loadDielectric('Diamond');
  // etc..
}

Scene.prototype.getMinScale = function() { return 1.0e-4; /* raymarch tolerance */ }
Scene.prototype.getMaxScale = function() { return 1.0e2; /* raymarch infinity */ }
Scene.prototype.envMap = function()  { return '<env map url>'; }
Scene.prototype.initGui = function(gui) { /* setup GUI */  }
Scene.prototype.syncShader = function(shader) { /* sync shader with GUI */ }
Scene.prototype.preframeCallback = function(snelly, gl) { /* custom logic */ }
Scene.prototype.postframeCallback = function(snelly, gl) { /* custom logic */ }

function onLoad() { snelly = new Snelly(new Scene()); animateLoop(); }
function animateLoop() { snelly.render(); window.requestAnimationFrame(animateLoop); }

</script>
</body>
```
The only mandatory function to implement in Scene is {@link Scene#shader}, the other are all optional. However the {@link Scene#init} function is almost always needed, to set the initial camera orientation at least.

### Geometry

A Snelly scene is assumed to consist of only (up to) three specified materials: a Metal, a Dielectric, and a plastic-like Surface ("uber" material). Each material has an associated surface which is defined by an SDF (signed distance function), i.e. where each function is negative corresponds to the interior of the body.

Thus we define the rendered scene geometry by specifying, via the {@link Scene#shader} call, three GLSL functions:

	- SDF_SURFACE(vec3 X): the SDF of the uber-surface material
	- SDF_METAL(vec3 X): the SDF of the (selected) metal material
	- SDF_DIELECTRIC(vec3 X): the SDF of the (selected) dielectric material
  
 We use [dat.GUI](https://workshop.chromeexperiments.com/examples/gui/#1--Basic-Usage) to provide a simple interactive UI for the scene and renderer state. Basic scene controls can be added via the 
 UI or animation control over the scene can be coded by adding uniform variables in the SDF functions, and setting them to the corresponding UI values in the {@link Scene#syncShader} function.

The details of the properties of the three material types can then be specified in {@link Scene#init} via the {@link Materials} object. Additional spatial dependence of the material surface properties can be introduced by providing modulating GLSL functions.

Procedural camera motion and scene animation can be authored (programmatically) via the pre- and post-frame callbacks.

As a standalone web page, a Snelly scene can be easily shared, for example by keeping the HTML file in a GitHub repository and simply linking to the file via [RawGit](https://rawgit.com/). 

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
	}
```

Other such env-maps are available from [here](https://github.com/portsmouth/envmaps) (convert to RawGit links first).

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
this.parameters.foo2 = ${this.parameters.foo2};
this.parameters.bar = ${this.parameters.bar};
this.frame = ${this.frame};
	`; 
}
```
With this code in place, the output on pressing 'O' is then a faithful representation of the entire scene state.



### Callbacks and animation

For implementation of custom animation logic, we use the simple mechanism of pre- and post-frame user callbacks, wherein the user can implement whatever logic he needs to programmatically animate the scene, camera, and materials. See the provided examples for details of how to use this implement animating scenes, and movie rendering.

```
Scene.prototype.preframeCallback = function(snelly, gl);
Scene.prototype.postframeCallback = function(snelly, gl);
```

