
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
    controls.keyPanSpeed = 1;
    controls.flySpeed = 0.01;

    // Renderer settings
    //renderer.width = 1280; // (if either width or height are not specified, render size will be taken from window
    //renderer.height = 720; // and will then auto-resize with the window)
    renderer.renderMode = 'pt';  // The other modes are: 'ao', 'normals'
    renderer.maxBounces = 9;
    renderer.maxMarchSteps = 512;
    renderer.radianceClamp = 0.4355179704016914; // (log scale)
    renderer.skyPower = 4.0;
    renderer.skyTintUp = [1.0, 1.0, 1.0];
    renderer.skyTintDown = [1.0, 1.0, 1.0];
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
* state to be tweaked in the UI then saved by copy-pasting code into the init function.
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
      return 'https://cdn.jsdelivr.net/gh/portsmouth/envmaps/74e9d389/HDR_040_Field_Bg.jpg';
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

* The code *must* define at least one of the four functions:
*```glsl
*      float SDF_SURFACE(vec3 X);
*      float SDF_METAL(vec3 X);
*      float SDF_DIELECTRIC(vec3 X);
*      float SDF_VOLUME(vec3 X);
*```
* Only the SDFs which are defined will be rendered.
*
* Optionally, any of the following functions defining the spatial dependence of material reflectances and roughnesses can be defined.
* The UI-exposed reflectance or roughness is supplied, and can be modified arbitrarily based on the supplied data of the primary ray hit point
* (and/or any other computed shader variables). The arguments to these functions are as follows:
*   - *refl_ui*: The UI-exposed constant reflectance color ([0,1] floats in sRGB color space)
*   - *roughness_ui*: The UI-exposed constant roughness
*   - *X*: world space hit point
*   - *N*: world space (outward) normal at the hit point
*   - *V*: world space view direction at the hit point (i.e. direction from the hit point to the eye).
*
* Note that the vec3 color returned is also in sRGB color space.
* (Any of these functions can be omitted, in which case they will be replaced with the default indicated).
*```glsl
        // return surface diffuse reflectance (defaults to just return the input UI constant refl_ui)
        vec3 SURFACE_DIFFUSE_REFLECTANCE(in vec3 refl_ui, in vec3 X, in vec3 N, in vec3 V);

        // return surface specular reflectance (defaults to just return the input UI constant refl_ui)
        vec3 SURFACE_SPECULAR_REFLECTANCE(in vec3 refl_ui, in vec3 X, in vec3 N, in vec3 V);

        // return surface roughness in [0,1] (defaults to just return the input roughness_ui)
        float SURFACE_ROUGHNESS(in float roughness_ui, in vec3 X, in vec3 N);

        // return local space normal (z is up)
        vec3 SURFACE_NORMAL_MAP(in vec3 X);

        // return metal roughness in [0,1] (defaults to just return the input UI constant roughness_ui)
        float METAL_ROUGHNESS(in float roughness_ui, in vec3 X, in vec3 N);

        // return metal specular reflectance (defaults to just return the input refl_ui)
        vec3 METAL_SPECULAR_REFLECTANCE(in vec3 refl_ui, in vec3 X, in vec3 N, in vec3 V);

        // return local space normal (z is up)
        vec3 METAL_NORMAL_MAP(in vec3 X);

        // return dielectric roughness in [0,1] (defaults to just return the input UI constant roughness_ui)
        float DIELECTRIC_ROUGHNESS(in float roughness_ui, in vec3 X, in vec3 N);

        // return dielectric specular reflectance (defaults to just return the input UI constant refl_ui)
        vec3 DIELECTRIC_SPECULAR_REFLECTANCE(in vec3 refl_ui, in vec3 X, in vec3 N, in vec3 V);

        // return local space normal (z is up)
        vec3 DIELECTRIC_NORMAL_MAP(in vec3 X);

        float VOLUME_EXTINCTION(float extinction_ui, in vec3 X);
        float VOLUME_EXTINCTION_MAX(float extinction_ui);
        vec3 VOLUME_SCATTERING_COLOR(in vec3 scattering_color_ui, in vec3 X);
        vec3 VOLUME_ABSORPTION_COLOR(in vec3 scattering_color_ui, in vec3 X);
        vec3 VOLUME_EMISSION(in vec3 emission_ui, in vec3 X);
*```
*Optionally, an init function can also be provided, which will be called first by each primary ray. 
*This is occasionally useful to prepare global variables for use during the succeeding computation for this pixel.
*```glsl
    void INIT();
*```
* @returns {String}
*/
Scene.prototype.shader = function()
{
    return ``;
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
* Optional callack which, if implemented, the renderer consults before
* each frame to determine whether to render. This is needed if the scene has to do some
* pre-processing, or load something, before rendering is possible.
* @param {Snelly} snelly - The snelly object
* @returns {boolean} - whether the scene is ready to render
*/
Scene.prototype.isReady = function(snelly)            
{
    return true;
}

/**
* Optional. Called whenever the UI is changed,
/* and must sync the params of the shader with the current UI settings
* @param {Snelly} snelly - The snelly object
* @param {GLU.this.Shader} shader - wrapper of webGL fragment shader, see {@link GLU.this.Shader}
*/
Scene.prototype.syncShader = function(snelly, shader)
{
    shader.uniformF("_foo", this.parameters.foo);
    shader.uniformF("_foo2", this.parameters.foo2);
    shader.uniformF("_bar", this.parameters.bar);
}


/**
* Optional. Gives the raytracer some indication of the (rough) typical length scale of this scene, 
* so it can set tolerances and defaults appropriately. 
* Defaults to 1.0.
* @returns {number}
*/
Scene.prototype.getLengthScale = function()
{
    return 1.0;
}

/**
* Optional. Gives the raytracer some indication of the (rough) minimum length scale, 
* so it can set tolerances appropriately. This sets the rough length scale of the smallest 
* resolvable structure. (Note that decreasing this will usually lead to longer render times).
* Defaults to 0.0001 of the scene length scale.
* @returns {number}
*/
Scene.prototype.getMinLengthScale = function()
{
    return 1.0e-4 * this.getLengthScale();
}


/**
* Optional. Gives the raytracer some indication of the (rough) maximum length scale, 
* so it can set tolerances appropriately. The raymarcher will march no further
* from the camera than this scale, thus it acts as the "far plane" distance.
* (Note that increasing this will usually lead to longer render times).
* Defaults to 100.0 of the scene length scale;
* @returns {number}
*/
Scene.prototype.getMaxLengthScale = function()
{
    return 1.0e2 * this.getLengthScale();
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
    return; 

    /*
    // The code here posts the framebuffer pixels to a local server, for sequence rendering.
    let renderer  = snelly.getRenderer();
    let targetSPP = 250.0;

    // User code to post webGL framebuffer data to local server for processing
    if (this.animFrame>=0 && renderer.getSPP()>=targetSPP && this.animFrame<=this.endFrame)
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
    */
}

/** 
 * Optional callback on key down.
 * @param {Event} Javascript keydown Event
 * @param {Snelly} snelly - The Snelly object
 * @param {WebGLRenderingContext} gl - The webGL context
 */
Scene.prototype.onkeydownCallback = function(event, snelly, gl)
{
    return;
}

