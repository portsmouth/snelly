
## Classes

<dl>
<dt><a href="#Snelly">Snelly</a></dt>
<dd></dd>
<dt><a href="#Scene">Scene</a></dt>
<dd></dd>
<dt><a href="#Renderer">Renderer</a></dt>
<dd></dd>
<dt><a href="#Material">Material</a></dt>
<dd></dd>
<dt><a href="#Surface">Surface</a> ⇐ <code><a href="#Material">Material</a></code></dt>
<dd></dd>
<dt><a href="#Volume">Volume</a> ⇐ <code><a href="#Material">Material</a></code></dt>
<dd></dd>
<dt><a href="#Metal">Metal</a> ⇐ <code><a href="#Material">Material</a></code></dt>
<dd></dd>
<dt><a href="#Dielectric">Dielectric</a> ⇐ <code><a href="#Material">Material</a></code></dt>
<dd></dd>
<dt><a href="#Materials">Materials</a></dt>
<dd></dd>
</dl>

## Objects

<dl>
<dt><a href="#GLU">GLU</a> : <code>object</code></dt>
<dd><p>Namespace for webGL utility wrappers.
Functions for loading shader uniform variables are exposed to the user
for convenience.</p>
</dd>
</dl>

<a name="Snelly"></a>

## Snelly
**Kind**: global class  

* [Snelly](#Snelly)
    * [new Snelly(sceneObj)](#new_Snelly_new)
    * [.getVersion()](#Snelly+getVersion) ⇒ <code>Array</code>
    * [.getRenderer()](#Snelly+getRenderer) ⇒ [<code>Renderer</code>](#Renderer)
    * [.getGUI()](#Snelly+getGUI) ⇒ <code>GUI</code>
    * [.getCamera()](#Snelly+getCamera) ⇒ <code>THREE.PerspectiveCamera</code>
    * [.getControls()](#Snelly+getControls) ⇒ <code>THREE.OrbitControls</code>
    * [.showGUI(show)](#Snelly+showGUI)
    * [.setStatus(statusText)](#Snelly+setStatus)
    * [.getGLContext()](#Snelly+getGLContext) ⇒ <code>WebGLRenderingContext</code>
    * [.getMaterials()](#Snelly+getMaterials) ⇒ [<code>Materials</code>](#Materials)
    * [.getSurface()](#Snelly+getSurface) ⇒ [<code>Surface</code>](#Surface)
    * [.getVolume()](#Snelly+getVolume) ⇒ [<code>Volume</code>](#Volume)
    * [.getUserTextureUnitStart()](#Snelly+getUserTextureUnitStart) ⇒ <code>number</code>

<a name="new_Snelly_new"></a>

### new Snelly(sceneObj)
Snelly is the global object providing access to all functionality in the system.


| Param | Type | Description |
| --- | --- | --- |
| sceneObj | [<code>Scene</code>](#Scene) | The user-defined scene |

<a name="Snelly+getVersion"></a>

### snelly.getVersion() ⇒ <code>Array</code>
Returns the current version number of the snelly system, in the format [1, 2, 3] (i.e. major, minor, patch version)

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getRenderer"></a>

### snelly.getRenderer() ⇒ [<code>Renderer</code>](#Renderer)
Access to the Renderer object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getGUI"></a>

### snelly.getGUI() ⇒ <code>GUI</code>
Access to the GUI object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getCamera"></a>

### snelly.getCamera() ⇒ <code>THREE.PerspectiveCamera</code>
Access to the camera object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
**Returns**: <code>THREE.PerspectiveCamera</code> - .  
<a name="Snelly+getControls"></a>

### snelly.getControls() ⇒ <code>THREE.OrbitControls</code>
Access to the camera controller object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+showGUI"></a>

### snelly.showGUI(show)
Programmatically show or hide the dat.GUI UI

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  

| Param | Type | Description |
| --- | --- | --- |
| show | <code>Boolean</code> | toggle |

<a name="Snelly+setStatus"></a>

### snelly.setStatus(statusText)
Specify arbitrary status text (one line only) to display in the lower right of the viewport

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  

| Param | Type | Description |
| --- | --- | --- |
| statusText | <code>Boolean</code> | text to display |

<a name="Snelly+getGLContext"></a>

### snelly.getGLContext() ⇒ <code>WebGLRenderingContext</code>
**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
**Returns**: <code>WebGLRenderingContext</code> - The webGL context  
<a name="Snelly+getMaterials"></a>

### snelly.getMaterials() ⇒ [<code>Materials</code>](#Materials)
Get materials object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getSurface"></a>

### snelly.getSurface() ⇒ [<code>Surface</code>](#Surface)
Get Surface object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getVolume"></a>

### snelly.getVolume() ⇒ [<code>Volume</code>](#Volume)
Get Volume object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getUserTextureUnitStart"></a>

### snelly.getUserTextureUnitStart() ⇒ <code>number</code>
**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
**Returns**: <code>number</code> - - the minimum texture unit for user supplied textures in the shader  
<a name="Scene"></a>

## Scene
**Kind**: global class  

* [Scene](#Scene)
    * [.init(snelly)](#Scene+init)
    * [.initGenerator()](#Scene+initGenerator)
    * [.envMap()](#Scene+envMap) ⇒ <code>String</code>
    * [.getName()](#Scene+getName) ⇒ <code>String</code>
    * [.getURL()](#Scene+getURL) ⇒ <code>String</code>
    * [.shader()](#Scene+shader) ⇒ <code>String</code>
    * [.initGui(gui)](#Scene+initGui)
    * [.isReady(snelly)](#Scene+isReady) ⇒ <code>boolean</code>
    * [.syncShader(snelly, shader)](#Scene+syncShader)
    * [.getLengthScale()](#Scene+getLengthScale) ⇒ <code>number</code>
    * [.getMinLengthScale()](#Scene+getMinLengthScale) ⇒ <code>number</code>
    * [.getMaxLengthScale()](#Scene+getMaxLengthScale) ⇒ <code>number</code>
    * [.preframeCallback(snelly, gl)](#Scene+preframeCallback)
    * [.postframeCallback(snelly, gl)](#Scene+postframeCallback)
    * [.onkeydownCallback(Javascript, snelly, gl)](#Scene+onkeydownCallback)

<a name="Scene+init"></a>

### scene.init(snelly)
Optionally (but usually), provide this function to set scene and renderer initial state.This is called only once during execution, on loading the scene HTML page (or on global reset via 'R' key).

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| snelly | [<code>Snelly</code>](#Snelly) | The snelly object |

<a name="Scene+initGenerator"></a>

### scene.initGenerator()
Optionally, provide this function which generates the init code to re-generate the current UI parameter settings. This will be dumped to the console (along with the rest of the UI state) on pressing key 'O', allowing the scene and rendererstate to be tweaked in the UI then saved by copy-pasting code into the init function.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+envMap"></a>

### scene.envMap() ⇒ <code>String</code>
Optionally, supply an env-map texture URL (must be a lat-long format image).(If this is function not implemented, or it returns the empty string, a uniformtemperature blackbody sky is used).

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getName"></a>

### scene.getName() ⇒ <code>String</code>
Optional name (displayed in UI)

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getURL"></a>

### scene.getURL() ⇒ <code>String</code>
Optional clickable URL (displayed in UI)

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+shader"></a>

### scene.shader() ⇒ <code>String</code>
Returns a chunk of GLSL code defining the SDFs which determine the geometry of uber-surface, metal and dielectric materials in the scene.Define also (optionally) functions giving the 3d spatial dependence of the material parameters.This function is mandatory!The code *must* define at least one of the four functions:```glsl     float SDF_SURFACE(vec3 X);     float SDF_METAL(vec3 X);     float SDF_DIELECTRIC(vec3 X);     float SDF_VOLUME(vec3 X);```Only the SDFs which are defined will be rendered.Optionally, any of the following functions defining the spatial dependence of material reflectances and roughnesses can be defined.The UI-exposed reflectance or roughness is supplied, and can be modified arbitrarily based on the supplied data of the primary ray hit point(and/or any other computed shader variables). The arguments to these functions are as follows:  - *refl_ui*: The UI-exposed constant reflectance color ([0,1] floats in sRGB color space)  - *roughness_ui*: The UI-exposed constant roughness  - *X*: world space hit point  - *N*: world space (outward) normal at the hit point  - *V*: world space view direction at the hit point (i.e. direction from the hit point to the eye).Note that the vec3 color returned is also in sRGB color space.(Any of these functions can be omitted, in which case they will be replaced with the default indicated).```glsl
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
        vec3 VOLUME_EMISSION(in vec3 emission_ui, in vec3 X);```Optionally, an init function can also be provided, which will be called first by each primary ray. This is occasionally useful to prepare global variables for use during the succeeding computation for this pixel.```glsl
    void INIT();```

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+initGui"></a>

### scene.initGui(gui)
Optional. Set up gui and callbacks for this scene

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| gui | <code>GUI</code> | wrapper for dat.GUI object |

<a name="Scene+isReady"></a>

### scene.isReady(snelly) ⇒ <code>boolean</code>
Optional callack which, if implemented, the renderer consults beforeeach frame to determine whether to render. This is needed if the scene has to do somepre-processing, or load something, before rendering is possible.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
**Returns**: <code>boolean</code> - - whether the scene is ready to render  

| Param | Type | Description |
| --- | --- | --- |
| snelly | [<code>Snelly</code>](#Snelly) | The snelly object |

<a name="Scene+syncShader"></a>

### scene.syncShader(snelly, shader)
Optional. Called whenever the UI is changed,
/* and must sync the params of the shader with the current UI settings

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| snelly | [<code>Snelly</code>](#Snelly) | The snelly object |
| shader | [<code>this.Shader</code>](#GLU.this.Shader) | wrapper of webGL fragment shader, see [this.Shader](#GLU.this.Shader) |

<a name="Scene+getLengthScale"></a>

### scene.getLengthScale() ⇒ <code>number</code>
Optional. Gives the raytracer some indication of the (rough) typical length scale of this scene, so it can set tolerances and defaults appropriately. Defaults to 1.0.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getMinLengthScale"></a>

### scene.getMinLengthScale() ⇒ <code>number</code>
Optional. Gives the raytracer some indication of the (rough) minimum length scale, so it can set tolerances appropriately. This sets the rough length scale of the smallest resolvable structure. (Note that decreasing this will usually lead to longer render times).Defaults to 0.0001 of the scene length scale.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getMaxLengthScale"></a>

### scene.getMaxLengthScale() ⇒ <code>number</code>
Optional. Gives the raytracer some indication of the (rough) maximum length scale, so it can set tolerances appropriately. The raymarcher will march no furtherfrom the camera than this scale, thus it acts as the "far plane" distance.(Note that increasing this will usually lead to longer render times).Defaults to 100.0 of the scene length scale;

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+preframeCallback"></a>

### scene.preframeCallback(snelly, gl)
Optional callback before every frame.Animation rendering logic can be implemented here by updating the scene programmatically according to the global time since init.

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| snelly | [<code>Snelly</code>](#Snelly) | The Snelly object |
| gl | <code>WebGLRenderingContext</code> | The webGL context |

<a name="Scene+postframeCallback"></a>

### scene.postframeCallback(snelly, gl)
Optional callback after every frame.Animation rendering logic can be implemented here by updating the scene programmatically according to the global time since init.

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| snelly | [<code>Snelly</code>](#Snelly) | The Snelly object |
| gl | <code>WebGLRenderingContext</code> | The webGL context |

<a name="Scene+onkeydownCallback"></a>

### scene.onkeydownCallback(Javascript, snelly, gl)
Optional callback on key down.

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| Javascript | <code>Event</code> | keydown Event |
| snelly | [<code>Snelly</code>](#Snelly) | The Snelly object |
| gl | <code>WebGLRenderingContext</code> | The webGL context |

<a name="Renderer"></a>

## Renderer
**Kind**: global class  
**Properties**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| width | <code>number</code> |  | (if not specified, fits to window) |
| height | <code>number</code> |  | (if not specified, fits to window) |
| [renderMode] | <code>String</code> | <code>&#x27;pt&#x27;</code> | rendering mode (either 'pt', 'ao', 'normals') |
| [dispersive] | <code>number</code> | <code>false</code> | enable dispersive (i.e. spectral) rendering |
| [maxSamplesPerFrame] | <code>number</code> | <code>1</code> | maximum number of per-pixel samples per frame |
| [maxSpp] | <code>number</code> | <code>1</code> | maximum number of samples-per-pixel, after which the render terminates |
| [maxBounces] | <code>number</code> | <code>3</code> | maximum number of surface bounces |
| [maxMarchSteps] | <code>number</code> | <code>256</code> | maximum number of raymarching steps per path segment |
| [maxStepsIsMiss] | <code>number</code> | <code>true</code> | whether rays which exceed max step count are considered hits or misses |
| [interactive] | <code>number</code> | <code>true</code> | if enabled, tries to maintain interactive frame rate at the expense of more noise |
| [goalFPS] | <code>number</code> | <code>10.0</code> | sampling will adjust to try to match goal FPS |
| [minsSPPToRedraw] | <code>number</code> | <code>0.0</code> | if >0.0, renderer will not redraw until the specified SPP have been accumulated |
| [radianceClamp] | <code>number</code> | <code>3.0</code> | clamp radiance to (10^) this max value, for firefly reduction |
| [wavelengthSamples] | <code>number</code> | <code>256</code> | number of samples to take over visible wavelength range |
| [exposure] | <code>number</code> | <code>0.0</code> | exposure, on a log scale |
| [gamma] | <code>number</code> | <code>2.2</code> | display gamma correction |
| [contrast] | <code>number</code> | <code>1.0</code> | tonemapping contrast |
| [saturation] | <code>number</code> | <code>1.0</code> | tonemapping saturation |
| [skyPower] | <code>number</code> | <code>4.0</code> | sky power (arbitrary units) |
| [skyTemperature] | <code>number</code> | <code>6000</code> | sky emission blackbody temperature (in Kelvin), used in dispersive mode only |
| [skyTintUp] | <code>Array</code> |  | sky color upwards tint |
| [skyTintDown] | <code>Array</code> |  | sky color downwards tint |
| [envMapVisible] | <code>number</code> | <code>true</code> | whether env map is visible to primary rays (otherwise black) |
| [envMapPhiRotation] | <code>number</code> | <code>0.0</code> | env map rotation about pole in degrees (0 to 360) |
| [envMapThetaRotation] | <code>number</code> | <code>0.0</code> | env map rotation about equator in degrees (0 to 180) |
| [envMapTransitionAngle] | <code>number</code> | <code>0.0</code> | angle over which env map tint transitions from upwards to downwards tint [degrees] |
| [sunPower] | <code>number</code> | <code>1.0</code> | sun power (arbitrary units) |
| [sunColor] | <code>Array</code> |  | sun color |
| [sunAngularSize] | <code>number</code> | <code>5.0</code> | sun angular size (degrees) |
| [sunLatitude] | <code>number</code> | <code>50.0</code> | sun latitude (degrees) |
| [sunLongitude] | <code>number</code> | <code>0.0</code> | sun longitude (degrees) |
| [sunVisibleDirectly] | <code>number</code> | <code>true</code> | whether sun is directly visible |
| [sphereLightRadius] | <code>number</code> | <code>0.0</code> | sphere light radius |
| [sphereLightPower] | <code>number</code> |  | sphere light power (arbitrary units) |
| [sphereLightPosition] | <code>Array</code> |  | whether sun is directly visible |
| [shadowStrength] | <code>number</code> | <code>1.0</code> | if <1.0, areas in shadow are not completely dark (provided mostly to allow rendering of occluded areas, e.g. fractals) |


* [Renderer](#Renderer)
    * [new Renderer()](#new_Renderer_new)
    * [.focusAt](#Renderer+focusAt)
        * [new Renderer.prototype.focusAt()](#new_Renderer+focusAt_new)
    * [.reset([no_recompile])](#Renderer+reset)
    * [.getSPP()](#Renderer+getSPP) ⇒ <code>number</code>

<a name="new_Renderer_new"></a>

### new Renderer()
Interface to the renderer. The rendering modes available are: - 'pt': pathtracer (uni-directional) - 'ao': ambient occlusion, colored via [Surface](#Surface) material diffuse albedo modulated by the `SURFACE_DIFFUSE_REFLECTANCE` shader function - 'normals': view normal at first hit as a color

<a name="Renderer+focusAt"></a>

### renderer.focusAt
**Kind**: instance class of [<code>Renderer</code>](#Renderer)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| x | <code>number</code> | x coordinate of focus as a fraction of current frame width |
| y | <code>number</code> | y coordinate of focus as a fraction of current frame width |

<a name="new_Renderer+focusAt_new"></a>

#### new Renderer.prototype.focusAt()
Adjust focal length to match the distance to the (first) object visible at a given location in the frame.

<a name="Renderer+reset"></a>

### renderer.reset([no_recompile])
Restart accumulating samples.

**Kind**: instance method of [<code>Renderer</code>](#Renderer)  

| Param | Type | Default | Description |
| --- | --- | --- | --- |
| [no_recompile] | <code>Boolean</code> | <code>false</code> | set to true if shaders need recompilation too |

<a name="Renderer+getSPP"></a>

### renderer.getSPP() ⇒ <code>number</code>
Read access to the current per-pixel sample count average.

**Kind**: instance method of [<code>Renderer</code>](#Renderer)  
**Returns**: <code>number</code> - - the current sample count average  
<a name="Material"></a>

## Material
**Kind**: global class  
<a name="new_Material_new"></a>

### new Material()
Generic material.

<a name="Surface"></a>

## Surface ⇐ [<code>Material</code>](#Material)
**Kind**: global class  
**Extends**: [<code>Material</code>](#Material)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| roughness | <code>number</code> | The surface roughness |
| ior | <code>number</code> | The surface coating ior |
| diffuseAlbedo | <code>Array</code> | The surface diffuse (RGB) color |
| specAlbedo | <code>Array</code> | The surface spec (RGB) color |

<a name="new_Surface_new"></a>

### new Surface()
Generic uber-surface material. Control via properties:

**Example**  
```js
surface.roughness = 0.05;surface.ior = 1.3530655391120507;surface.diffuseAlbedo = [0.5, 0.5, 0.5];surface.specAlbedo = [0.0, 0.0, 0.0];
```
<a name="Volume"></a>

## Volume ⇐ [<code>Material</code>](#Material)
**Kind**: global class  
**Extends**: [<code>Material</code>](#Material)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| mfp | <code>number</code> | MFP in units of inverse scene scale (gives grey extinction as inverse MFP) |
| maxOpticalDepth | <code>number</code> | maximum optical depth (in any channel), used to bound attenuation to infinity |
| scatteringColor | <code>Array</code> | Scattering (RGB) color (multiplies grey extinction to give per-channel scattering coefficient) |
| absorptionColor | <code>Array</code> | The absorption (RGB) color (multiplies extinction to give per-channel absorption coefficient) |
| anisotropy | <code>number</code> | Phase function anisotropy in [-1,1] |
| emission | <code>number</code> | emission power magnitude |
| emissionColor | <code>Array</code> | emission color (multiplies emission to give per-channel emission) |

<a name="new_Volume_new"></a>

### new Volume()
Volumetric material (a homogeneous atmosphere, with heterogeneous emission). Control via properties:

**Example**  
```js
volume.mfp = 0.1;volume.scatteringColor = [0.5, 0.5, 0.5];volume.absorptionColor = [0.0, 0.5, 0.0];volume.anisotropy = 0.0;volume.emission = 0.0;volume.emissionColor = [0.5, 0.5, 0.5];
```
<a name="Metal"></a>

## Metal ⇐ [<code>Material</code>](#Material)
**Kind**: global class  
**Extends**: [<code>Material</code>](#Material)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| roughness | <code>number</code> | The metal surface roughness |

<a name="new_Metal_new"></a>

### new Metal()
Generic metal material. Supported physical metals are:``` "Aluminium" "Brass" "Calcium" "Chromium" "Cobalt" "Copper"  "Gold"    "Iridium" "Iron"  "Lead"    "Mercury" "Molybdenum" "Nickel" "Palladium" "Platinum" "Silicon" "Silver" "Titanium" "Tungsten" "Vanadium" "Zinc" "Zirconium"```

**Example**  
```js
let metal = materials.loadMetal('Gold');metal.roughness = 0.05;
```
<a name="Dielectric"></a>

## Dielectric ⇐ [<code>Material</code>](#Material)
**Kind**: global class  
**Extends**: [<code>Material</code>](#Material)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| roughness | <code>number</code> | The dielectric surface roughness |
| absorptionColor | <code>array</code> | The dielectric surface absorption color |
| absorptionScale | <code>number</code> | The dielectric surface absorption scale (m.f.p in multiples of scene max scale) |

<a name="new_Dielectric_new"></a>

### new Dielectric()
Generic dielectric material. Supported physical dielectrics are:```glsl "Constant IOR dielectric" "Glass (BK7)" "Glass (K7)" "Glass (F5)" "Glass (LAFN7)" "Glass (LASF35)" "Glass (N-LAK33A)" "Glass (N-FK51A)" "Glass (SF4)" "Glass (SF67)" "Water" "Polycarbonate" "Glycerol" "Liquid Crystal (E7)" "Diamond" "Quartz" "Fused Silica" "Sapphire" "Sodium Chloride" "Proustite" "Rutile" "Silver Chloride"```

**Example**  
```js
let dielectric = materials.loadDielectric('Diamond');dielectric.absorptionColor = [1.0, 1.0, 1.0];dielectric.absorptionScale = 1.0; // mfp in multiples of scene scaledielectric.roughness = 0.030443974630021145;
```
<a name="Materials"></a>

## Materials
**Kind**: global class  

* [Materials](#Materials)
    * [new Materials()](#new_Materials_new)
    * [.loadDielectric(dielectricName)](#Materials+loadDielectric) ⇒ [<code>Dielectric</code>](#Dielectric)
    * [.getDielectric()](#Materials+getDielectric) ⇒ [<code>Dielectric</code>](#Dielectric)
    * [.loadMetal(metalName)](#Materials+loadMetal) ⇒ [<code>Metal</code>](#Metal)
    * [.getMetal()](#Materials+getMetal) ⇒ [<code>Metal</code>](#Metal)
    * [.getSurface()](#Materials+getSurface) ⇒ [<code>Surface</code>](#Surface)
    * [.getVolume()](#Materials+getVolume) ⇒ [<code>Surface</code>](#Surface)

<a name="new_Materials_new"></a>

### new Materials()
This object controls the properties of the three basic material types: - Dielectric (multiple different sub-types) - Metal (multiple different sub-types) - Surface (an uber-shader like materal)

<a name="Materials+loadDielectric"></a>

### materials.loadDielectric(dielectricName) ⇒ [<code>Dielectric</code>](#Dielectric)
Load the desired Dielectric object by name. Supported dielectrics are:```glsl "Constant IOR dielectric" "Glass (BK7)" "Glass (K7)" "Glass (F5)" "Glass (LAFN7)" "Glass (LASF35)" "Glass (N-LAK33A)" "Glass (N-FK51A)" "Glass (SF4)" "Glass (SF67)" "Water" "Polycarbonate" "Glycerol" "Liquid Crystal (E7)" "Diamond" "Quartz" "Fused Silica" "Sapphire" "Sodium Chloride" "Proustite" "Rutile" "Silver Chloride"```

**Kind**: instance method of [<code>Materials</code>](#Materials)  
**Returns**: [<code>Dielectric</code>](#Dielectric) - - the loaded dielectric  

| Param | Type | Description |
| --- | --- | --- |
| dielectricName | <code>String</code> | one of the names listed above |

<a name="Materials+getDielectric"></a>

### materials.getDielectric() ⇒ [<code>Dielectric</code>](#Dielectric)
Get the currently loaded Dielectric object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="Materials+loadMetal"></a>

### materials.loadMetal(metalName) ⇒ [<code>Metal</code>](#Metal)
Load the desired Metal object by name. Supported metals are:``` "Aluminium" "Brass" "Calcium" "Chromium" "Cobalt"   "Copper"   "Gold"     "Iridium" "Iron"     "Lead"     "Mercury"  "Molybdenum" "Nickel" "Palladium" "Platinum" "Silicon" "Silver"  "Titanium" "Tungsten" "Vanadium" "Zinc"  "Zirconium"```

**Kind**: instance method of [<code>Materials</code>](#Materials)  
**Returns**: [<code>Metal</code>](#Metal) - - the loaded metal  

| Param | Type | Description |
| --- | --- | --- |
| metalName | <code>String</code> | one of the names listed above |

<a name="Materials+getMetal"></a>

### materials.getMetal() ⇒ [<code>Metal</code>](#Metal)
Get the currently loaded Metal object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="Materials+getSurface"></a>

### materials.getSurface() ⇒ [<code>Surface</code>](#Surface)
Get the Surface object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="Materials+getVolume"></a>

### materials.getVolume() ⇒ [<code>Surface</code>](#Surface)
Get the Volume object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="GLU"></a>

## GLU : <code>object</code>
Namespace for webGL utility wrappers.Functions for loading shader uniform variables are exposed to the userfor convenience.

**Kind**: global namespace  

* [GLU](#GLU) : <code>object</code>
    * [.this.Shader](#GLU.this.Shader)
        * [new this.Shader()](#new_GLU.this.Shader_new)
        * [.uniformI(name, i)](#GLU.this.Shader.uniformI)
        * [.uniformF(name, f)](#GLU.this.Shader.uniformF)
        * [.uniform2F(name, f1, f2)](#GLU.this.Shader.uniform2F)
        * [.uniform1Fv(name, fvec)](#GLU.this.Shader.uniform1Fv)
        * [.uniform2Fv(name, fvec2)](#GLU.this.Shader.uniform2Fv)
        * [.uniform3F(name, f1, f2, f3)](#GLU.this.Shader.uniform3F)
        * [.uniform3Fv(name, fvec3)](#GLU.this.Shader.uniform3Fv)
        * [.uniform4F(name, f1, f2, f3, f4)](#GLU.this.Shader.uniform4F)
        * [.uniform4Fv(name, fvec4)](#GLU.this.Shader.uniform4Fv)
        * [.uniformMatrix4fv(name, matrixArray16)](#GLU.this.Shader.uniformMatrix4fv)

<a name="GLU.this.Shader"></a>

### GLU.this.Shader
**Kind**: static class of [<code>GLU</code>](#GLU)  

* [.this.Shader](#GLU.this.Shader)
    * [new this.Shader()](#new_GLU.this.Shader_new)
    * [.uniformI(name, i)](#GLU.this.Shader.uniformI)
    * [.uniformF(name, f)](#GLU.this.Shader.uniformF)
    * [.uniform2F(name, f1, f2)](#GLU.this.Shader.uniform2F)
    * [.uniform1Fv(name, fvec)](#GLU.this.Shader.uniform1Fv)
    * [.uniform2Fv(name, fvec2)](#GLU.this.Shader.uniform2Fv)
    * [.uniform3F(name, f1, f2, f3)](#GLU.this.Shader.uniform3F)
    * [.uniform3Fv(name, fvec3)](#GLU.this.Shader.uniform3Fv)
    * [.uniform4F(name, f1, f2, f3, f4)](#GLU.this.Shader.uniform4F)
    * [.uniform4Fv(name, fvec4)](#GLU.this.Shader.uniform4Fv)
    * [.uniformMatrix4fv(name, matrixArray16)](#GLU.this.Shader.uniformMatrix4fv)

<a name="new_GLU.this.Shader_new"></a>

#### new this.Shader()
Represents a webGL vertex or fragment shader:

<a name="GLU.this.Shader.uniformI"></a>

#### this.Shader.uniformI(name, i)
Provide an integer (via uniform1i) to the currently bound shader

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| i | <code>number</code> | The integer value |

<a name="GLU.this.Shader.uniformF"></a>

#### this.Shader.uniformF(name, f)
Provide a float (via uniform1f) to the currently bound shader

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f | <code>number</code> | The float value |

<a name="GLU.this.Shader.uniform2F"></a>

#### this.Shader.uniform2F(name, f1, f2)
Provide a vec2 uniform (via uniform2f) to the currently bound shader

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f1 | <code>number</code> | The first float value |
| f2 | <code>number</code> | The second float value |

<a name="GLU.this.Shader.uniform1Fv"></a>

#### this.Shader.uniform1Fv(name, fvec)
Provide an array of floats (via uniform1Fv) to the currently bound shader  i.e. the shader declares e.g. `uniform float values[19];`

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec | <code>Float32Array</code> | An array of floats |

<a name="GLU.this.Shader.uniform2Fv"></a>

#### this.Shader.uniform2Fv(name, fvec2)
Provide an array of vec2 (via uniform2fv) to the currently bound shader  i.e. the shader declares e.g. `uniform vec2 vectors[19];`

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec2 | <code>Float32Array</code> | An array of floats, 2 per vector |

<a name="GLU.this.Shader.uniform3F"></a>

#### this.Shader.uniform3F(name, f1, f2, f3)
Provide a vec3 uniform (via uniform3f) to the currently bound shader

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f1 | <code>number</code> | The first float value |
| f2 | <code>number</code> | The second float value |
| f3 | <code>number</code> | The third float value |

<a name="GLU.this.Shader.uniform3Fv"></a>

#### this.Shader.uniform3Fv(name, fvec3)
Provide an array of vec3 (via uniform3fv) to the currently bound shader  i.e. the shader declares e.g. `uniform vec3 vectors[19];`

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec3 | <code>Float32Array</code> | An array of floats, 3 per vector |

<a name="GLU.this.Shader.uniform4F"></a>

#### this.Shader.uniform4F(name, f1, f2, f3, f4)
Provide a vec4 uniform (via uniform4F) to the currently bound shader

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f1 | <code>number</code> | The first float value |
| f2 | <code>number</code> | The second float value |
| f3 | <code>number</code> | The third float value |
| f4 | <code>number</code> | The fourth float value |

<a name="GLU.this.Shader.uniform4Fv"></a>

#### this.Shader.uniform4Fv(name, fvec4)
Provide an array of vec4 (via uniform4fv) to the currently bound shader  i.e. the shader declares e.g. `uniform vec4 vectors[19];`

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec4 | <code>Float32Array</code> | An array of floats, 4 per vector |

<a name="GLU.this.Shader.uniformMatrix4fv"></a>

#### this.Shader.uniformMatrix4fv(name, matrixArray16)
Provide a matrix (via uniformMatrix4fv) to the currently bound shader i.e. the shader declares e.g. `uniform mat4 matrix;`

**Kind**: static method of [<code>this.Shader</code>](#GLU.this.Shader)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| matrixArray16 | <code>Float32Array</code> | An array of 16 floats |


