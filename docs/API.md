
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
    * [.showGUI(showGUI)](#Snelly+showGUI)
    * [.getMaterials()](#Snelly+getMaterials) ⇒ [<code>Materials</code>](#Materials)
    * [.getSurface()](#Snelly+getSurface) ⇒ [<code>Surface</code>](#Surface)

<a name="new_Snelly_new"></a>

### new Snelly(sceneObj)
Snelly is the global object providing access to all functionality in the system.


| Param | Type | Description |
| --- | --- | --- |
| sceneObj | [<code>Scene</code>](#Scene) | The user-defined scene |

<a name="Snelly+getVersion"></a>

### snelly.getVersion() ⇒ <code>Array</code>
Returns the current version number of the snelly system, in the format [1, 2, 3]

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
**Returns**: <code>THREE.OrbitControls</code> - .  
<a name="Snelly+showGUI"></a>

### snelly.showGUI(showGUI)
Programmatically show or hide the dat.GUI UI

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  

| Param | Type | Description |
| --- | --- | --- |
| showGUI | <code>Boolean</code> | toggle |

<a name="Snelly+getMaterials"></a>

### snelly.getMaterials() ⇒ [<code>Materials</code>](#Materials)
Get materials object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
<a name="Snelly+getSurface"></a>

### snelly.getSurface() ⇒ [<code>Surface</code>](#Surface)
Get Surface object

**Kind**: instance method of [<code>Snelly</code>](#Snelly)  
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
    * [.initGui(The)](#Scene+initGui)
    * [.syncShader(The)](#Scene+syncShader)
    * [.getMinScale()](#Scene+getMinScale) ⇒ <code>number</code>
    * [.getMaxScale()](#Scene+getMaxScale) ⇒ <code>number</code>
    * [.preframeCallback(The, The)](#Scene+preframeCallback)
    * [.postframeCallback(The, The)](#Scene+postframeCallback)

<a name="Scene+init"></a>

### scene.init(snelly)
Optionally (but usually), provide this function to set scene and renderer initial state.

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| snelly | [<code>Snelly</code>](#Snelly) | The snelly object |

<a name="Scene+initGenerator"></a>

### scene.initGenerator()
Optionally, provide this function which generates the init code to re-generate 
the current UI parameter settings. This will be dumped to the console (along with 
the rest of the UI state) on pressing key 'O', allowing the scene and renderer
state to be tweaked in the UI then saved by copy-pasting code into the init function below.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+envMap"></a>

### scene.envMap() ⇒ <code>String</code>
Optionally, supply an env-map texture URL (must be a lat-long format image).
(If this is function not implemented, or it returns the empty string, a uniform
temperature blackbody sky is used).

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
Returns a chunk of GLSL code defining the SDFs which determine the geometry of uber-surface, metal and dielectric materials in the scene.
Define also (optionally) functions giving the 3d spatial dependence of the material parameters.
This function is mandatory!
The code *must* define at least one of the three functions:
```glsl
     float SDF_METAL(vec3 X);
     float SDF_DIELECTRIC(vec3 X);
     float SDF_DIFFUSE(vec3 X);
```
(If only one or two of these are present, the others will not be rendered).
Optionally, any of the following functions defining the spatial *modulation* of material parameters can be defined.
(Any of the user-defined functions below can be omitted, in which case they will be replaced with the default indicated).
```glsl

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
```

**Kind**: instance method of [<code>Scene</code>](#Scene)  
**Returns**: <code>String</code> - .  
<a name="Scene+initGui"></a>

### scene.initGui(The)
Optional. Set up gui and callbacks for this scene

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| The | <code>GUI</code> | GUI object |

<a name="Scene+syncShader"></a>

### scene.syncShader(The)
Optional. Called whenever the UI is changed,
/* and must sync the params of the shader with the current UI settings

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| The | <code>Shader</code> | Shader object |

<a name="Scene+getMinScale"></a>

### scene.getMinScale() ⇒ <code>number</code>
Optional. Gives the raytracer some indication of the (rough) minimum length scale, 
so it can set tolerances appropriately. This sets the rough length scale of the smallest 
resolvable structure. (Note that decreasing this will usually lead to longer render times).
Defaults to 0.0001.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getMaxScale"></a>

### scene.getMaxScale() ⇒ <code>number</code>
Optional. Gives the raytracer some indication of the (rough) maximum length scale, 
so it can set tolerances appropriately. The raymarcher will march no further
from the camera than this scale, thus it acts as the "far plane" distance.
(Note that increasing this will usually lead to longer render times).
Defaults to 100.0.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+preframeCallback"></a>

### scene.preframeCallback(The, The)
Optional callback before every frame.
Animation rendering logic can be implemented here by updating the scene 
programmatically according to the global time since init

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| The | [<code>Snelly</code>](#Snelly) | snelly object |
| The | <code>WebGLRenderingContext</code> | webGL context |

<a name="Scene+postframeCallback"></a>

### scene.postframeCallback(The, The)
Optional callback after every frame.
Animation rendering logic can be implemented here by updating the scene 
programmatically according to the global time since init

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| The | [<code>Snelly</code>](#Snelly) | snelly object |
| The | <code>WebGLRenderingContext</code> | webGL context |

<a name="Renderer"></a>

## Renderer
**Kind**: global class  
**Properties**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| width | <code>number</code> |  | (if not specified, fits to window) |
| height | <code>number</code> |  | (if not specified, fits to window) |
| renderMode | <code>String</code> | <code>&#x27;pt&#x27;</code> | rendering mode (either 'pt', 'ao', or 'normals') |
| maxMarchSteps | <code>number</code> | <code>256</code> | maximum number of raymarching steps per path segment |
| radianceClamp | <code>number</code> | <code>3.0</code> | clamp radiance to (10^) this max value, for firefly reduction |
| skyPower | <code>number</code> | <code>4.0</code> | sky power (arbitrary units) |
| skyTemperature | <code>number</code> | <code>6000</code> | sky temperature (in Kelvin) |
| exposure | <code>number</code> | <code>4.5</code> | exposure, on a log scale |
| gamma | <code>number</code> | <code>2.2</code> | display gamma correction |
| whitepoint | <code>number</code> | <code>2.0</code> | tonemapping whitepoint |
| goalFPS | <code>number</code> | <code>10.0</code> | sampling will adjust to try to match goal FPS |
| minsSPPToRedraw | <code>number</code> | <code>0.0</code> | if >0.0, renderer will not redraw until the specified SPP have been accumulated |


* [Renderer](#Renderer)
    * [new Renderer()](#new_Renderer_new)
    * [.reset([no_recompile])](#Renderer+reset)

<a name="new_Renderer_new"></a>

### new Renderer()
Interface to the pathtracer. The exposed properties and their defaults are:

<a name="Renderer+reset"></a>

### renderer.reset([no_recompile])
Restart accumulating samples.

**Kind**: instance method of [<code>Renderer</code>](#Renderer)  

| Param | Type | Default | Description |
| --- | --- | --- | --- |
| [no_recompile] | <code>Boolean</code> | <code>false</code> | set to true if shaders need recompilation too |

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
surface.roughness = 0.05;
surface.ior = 1.3530655391120507;
surface.diffuseAlbedo = [0.5, 0.5, 0.5];
surface.specAlbedo = [0.0, 0.0, 0.0];
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
Generic metal material. Supported physical metals are:
```
 "Aluminium"
 "Brass",   
 "Calcium", 
 "Chromium",
 "Cobalt",  
 "Copper",  
 "Gold",    
 "Iridium", 
 "Iron",    
 "Lead",    
 "Mercury", 
 "Molybdenum
 "Nickel",  
 "Palladium"
 "Platinum",
 "Silicon", 
 "Silver",  
 "Titanium",
 "Tungsten",
 "Vanadium",
 "Zinc",    
 "Zirconium"
```

**Example**  
```js
let metal = materials.loadMetal('Gold');
metal.roughness = 0.05;
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
Generic dielectric material. Supported physical dielectrics are:
```glsl
 "Constant IOR dielectric"
 "Glass (BK7)"
 "Glass (K7)"
 "Glass (F5)"
 "Glass (LAFN7)"
 "Glass (LASF35)"
 "Glass (N-LAK33A)"
 "Glass (N-FK51A)"
 "Glass (SF4)"
 "Glass (SF67)"
 "Water"
 "Polycarbonate"
 "Glycerol"
 "Liquid Crystal (E7)"
 "Diamond"
 "Quartz"
 "Fused Silica"
 "Sapphire"
 "Sodium Chloride"
 "Proustite"
 "Rutile"
 "Silver Chloride"
```

**Example**  
```js
let dielectric = materials.loadDielectric('Diamond');
dielectric.absorptionColor = [1.0, 1.0, 1.0];
dielectric.absorptionScale = 1.0; // mfp in multiples of scene scale
dielectric.roughness = 0.030443974630021145;
```
<a name="Materials"></a>

## Materials
**Kind**: global class  

* [Materials](#Materials)
    * [new Materials()](#new_Materials_new)
    * [.loadDielectric(dielectricName)](#Materials+loadDielectric) ⇒ [<code>Dielectric</code>](#Dielectric)
    * [.loadMetal(metalName)](#Materials+loadMetal) ⇒ [<code>Metal</code>](#Metal)
    * [.getDielectric()](#Materials+getDielectric) ⇒ [<code>Dielectric</code>](#Dielectric)
    * [.getMetal()](#Materials+getMetal) ⇒ [<code>Metal</code>](#Metal)
    * [.getSurface()](#Materials+getSurface) ⇒ [<code>Surface</code>](#Surface)

<a name="new_Materials_new"></a>

### new Materials()
This object controls the properties of the three basic material types:
 - Dielectric (multiple different sub-types)
 - Metal (multiple different sub-types)
 - Surface (an uber-shader like materal)

<a name="Materials+loadDielectric"></a>

### materials.loadDielectric(dielectricName) ⇒ [<code>Dielectric</code>](#Dielectric)
Load the desired Dielectric object by name. Supported dielectrics are:
```glsl
 "Constant IOR dielectric"
 "Glass (BK7)"
 "Glass (K7)"
 "Glass (F5)"
 "Glass (LAFN7)"
 "Glass (LASF35)"
 "Glass (N-LAK33A)"
 "Glass (N-FK51A)"
 "Glass (SF4)"
 "Glass (SF67)"
 "Water"
 "Polycarbonate"
 "Glycerol"
 "Liquid Crystal (E7)"
 "Diamond"
 "Quartz"
 "Fused Silica"
 "Sapphire"
 "Sodium Chloride"
 "Proustite"
 "Rutile"
 "Silver Chloride"
```

**Kind**: instance method of [<code>Materials</code>](#Materials)  
**Returns**: [<code>Dielectric</code>](#Dielectric) - - the loaded dielectric  

| Param | Type | Description |
| --- | --- | --- |
| dielectricName | <code>String</code> | one of the names listed above |

<a name="Materials+loadMetal"></a>

### materials.loadMetal(metalName) ⇒ [<code>Metal</code>](#Metal)
Load the desired Metal object by name. Supported metals are:
```
 "Aluminium"
 "Brass",   
 "Calcium", 
 "Chromium",
 "Cobalt",  
 "Copper",  
 "Gold",    
 "Iridium", 
 "Iron",    
 "Lead",    
 "Mercury", 
 "Molybdenum
 "Nickel",  
 "Palladium"
 "Platinum",
 "Silicon", 
 "Silver",  
 "Titanium",
 "Tungsten",
 "Vanadium",
 "Zinc",    
 "Zirconium"
```

**Kind**: instance method of [<code>Materials</code>](#Materials)  
**Returns**: [<code>Metal</code>](#Metal) - - the loaded metal  

| Param | Type | Description |
| --- | --- | --- |
| metalName | <code>String</code> | one of the names listed above |

<a name="Materials+getDielectric"></a>

### materials.getDielectric() ⇒ [<code>Dielectric</code>](#Dielectric)
Get the currently loaded Dielectric object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="Materials+getMetal"></a>

### materials.getMetal() ⇒ [<code>Metal</code>](#Metal)
Get the currently loaded Metal object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="Materials+getSurface"></a>

### materials.getSurface() ⇒ [<code>Surface</code>](#Surface)
Get the Surface object.

**Kind**: instance method of [<code>Materials</code>](#Materials)  
<a name="GLU"></a>

## GLU : <code>object</code>
Namespace for webGL utility wrappers.
Functions for loading shader uniform variables are exposed to the user
for convenience.

**Kind**: global namespace  

* [GLU](#GLU) : <code>object</code>
    * [.this.Shader](#GLU.this.Shader)
        * [new this.Shader()](#new_GLU.this.Shader_new)

<a name="GLU.this.Shader"></a>

### GLU.this.Shader
**Kind**: static class of [<code>GLU</code>](#GLU)  
<a name="new_GLU.this.Shader_new"></a>

#### new this.Shader()
Represents a webGL vertex or fragment shader:


