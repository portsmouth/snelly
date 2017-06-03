
## Classes

<dl>
<dt><a href="#Snelly is the global object providing access to all functionality in the system.">Snelly is the global object providing access to all functionality in the system.</a></dt>
<dd></dd>
<dt><a href="#Scene">Scene</a></dt>
<dd></dd>
<dt><a href="#Generic material.">Generic material.</a></dt>
<dd></dd>
<dt><a href="#Surface">Surface</a> ⇐ <code>Material</code></dt>
<dd></dd>
<dt><a href="#Metal">Metal</a> ⇐ <code>Material</code></dt>
<dd></dd>
<dt><a href="#Dielectric">Dielectric</a> ⇐ <code>Material</code></dt>
<dd></dd>
<dt><a href="#This object controls the properties of the three basic material types_
 - Dielectric (multiple different sub-types)
 - Metal (multiple different sub-types)
 - Surface (an uber-shader like materal)">This object controls the properties of the three basic material types:
 - Dielectric (multiple different sub-types)
 - Metal (multiple different sub-types)
 - Surface (an uber-shader like materal)</a></dt>
<dd></dd>
</dl>

## Objects

<dl>
<dt><a href="#GLU">GLU</a> : <code>object</code></dt>
<dd><p>Namespace for webGL utility wrappers.
Functions for loading uniform variables is exposed to the user
for convenience.</p>
</dd>
</dl>

<a name="Snelly is the global object providing access to all functionality in the system."></a>

## Snelly is the global object providing access to all functionality in the system.
**Kind**: global class  
<a name="new_Snelly is the global object providing access to all functionality in the system._new"></a>

### new Snelly is the global object providing access to all functionality in the system.(sceneObj)

| Param | Type | Description |
| --- | --- | --- |
| sceneObj | [<code>Scene</code>](#Scene) | The user-defined scene |

<a name="Scene"></a>

## Scene
**Kind**: global class  

* [Scene](#Scene)
    * [.init(snelly)](#Scene+init)
    * [.initGenerator()](#Scene+initGenerator)
    * [.envMap()](#Scene+envMap)
    * [.getName()](#Scene+getName) ⇒ <code>String</code>
    * [.getURL()](#Scene+getURL)
    * [.shader()](#Scene+shader) ⇒ <code>String</code>
    * [.initGui(The)](#Scene+initGui)
    * [.syncShader(The)](#Scene+syncShader)
    * [.getMinScale()](#Scene+getMinScale)
    * [.getMaxScale()](#Scene+getMaxScale)
    * [.preframeCallback(The, The)](#Scene+preframeCallback)
    * [.postframeCallback(The, The)](#Scene+postframeCallback)

<a name="Scene+init"></a>

### scene.init(snelly)
Optionally (but usually), provide this function to set scene and renderer initial state.

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| snelly | <code>Snelly</code> | The snelly object |

<a name="Scene+initGenerator"></a>

### scene.initGenerator()
Optionally, provide this function which generates the init code to re-generate 
the current UI parameter settings. This will be dumped to the console (along with 
the rest of the UI state) on pressing key 'O', allowing the scene and renderer
state to be tweaked in the UI then saved by copy-pasting code into the init function below.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+envMap"></a>

### scene.envMap()
Optionally, supply an env-map texture URL (must be a lat-long format image).
(If this is function not implemented, or it returns the empty string, a uniform
temperature blackbody sky is used).

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getName"></a>

### scene.getName() ⇒ <code>String</code>
Optional name (displayed in UI)

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getURL"></a>

### scene.getURL()
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

### scene.getMinScale()
Optional. Gives the raytracer some indication of the (rough) minimum length scale, 
so it can set tolerances appropriately. This sets the rough length scale of the smallest 
resolvable structure. (Note that decreasing this will usually lead to longer render times).
Defaults to 0.0001.

**Kind**: instance method of [<code>Scene</code>](#Scene)  
<a name="Scene+getMaxScale"></a>

### scene.getMaxScale()
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
| The | <code>Snelly</code> | snelly object |
| The | <code>WebGLRenderingContext</code> | webGL context |

<a name="Scene+postframeCallback"></a>

### scene.postframeCallback(The, The)
Optional callback after every frame.
Animation rendering logic can be implemented here by updating the scene 
programmatically according to the global time since init

**Kind**: instance method of [<code>Scene</code>](#Scene)  

| Param | Type | Description |
| --- | --- | --- |
| The | <code>Snelly</code> | snelly object |
| The | <code>WebGLRenderingContext</code> | webGL context |

<a name="Generic material."></a>

## Generic material.
**Kind**: global class  
<a name="Surface"></a>

## Surface ⇐ <code>Material</code>
**Kind**: global class  
**Extends**: <code>Material</code>  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| roughness | <code>number</code> | The surface roughness |
| ior | <code>number</code> | The surface coating ior |
| diffuseAlbedo | <code>Array</code> | The surface diffuse (RGB) color |
| specAlbedo | <code>Array</code> | The surface spec (RGB) color |

<a name="Metal"></a>

## Metal ⇐ <code>Material</code>
**Kind**: global class  
**Extends**: <code>Material</code>  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| roughness | <code>number</code> | The metal surface roughness |

<a name="Dielectric"></a>

## Dielectric ⇐ <code>Material</code>
**Kind**: global class  
**Extends**: <code>Material</code>  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| roughness | <code>number</code> | The dielectric surface roughness |

<a name="This object controls the properties of the three basic material types_
 - Dielectric (multiple different sub-types)
 - Metal (multiple different sub-types)
 - Surface (an uber-shader like materal)"></a>

## This object controls the properties of the three basic material types:
 - Dielectric (multiple different sub-types)
 - Metal (multiple different sub-types)
 - Surface (an uber-shader like materal)
**Kind**: global class  
<a name="GLU"></a>

## GLU : <code>object</code>
Namespace for webGL utility wrappers.
Functions for loading uniform variables is exposed to the user
for convenience.

**Kind**: global namespace  

* [GLU](#GLU) : <code>object</code>
    * [.this.Shader#uniformI(name, i)](#GLU.this.Shader+uniformI)
    * [.this.Shader#uniformF(name, f)](#GLU.this.Shader+uniformF)
    * [.this.Shader#uniform2F(name, f1, f2)](#GLU.this.Shader+uniform2F)
    * [.this.Shader#uniform1Fv(name, fvec)](#GLU.this.Shader+uniform1Fv)
    * [.this.Shader#uniform2Fv(name, fvec2)](#GLU.this.Shader+uniform2Fv)
    * [.this.Shader#uniform3F(name, f1, f2, f3)](#GLU.this.Shader+uniform3F)
    * [.this.Shader#uniform3Fv(name, fvec3)](#GLU.this.Shader+uniform3Fv)
    * [.this.Shader#uniform4F(name, f1, f2, f3, f4)](#GLU.this.Shader+uniform4F)
    * [.this.Shader#uniform4Fv(name, fvec4)](#GLU.this.Shader+uniform4Fv)
    * [.this.Shader#uniformMatrix4fv(name, matrixArray16)](#GLU.this.Shader+uniformMatrix4fv)

<a name="GLU.this.Shader+uniformI"></a>

### GLU.this.Shader#uniformI(name, i)
Provide an integer (via uniform1i) to the currently bound shader

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| i | <code>number</code> | The integer value |

<a name="GLU.this.Shader+uniformF"></a>

### GLU.this.Shader#uniformF(name, f)
Provide a float (via uniform1f) to the currently bound shader

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f | <code>number</code> | The float value |

<a name="GLU.this.Shader+uniform2F"></a>

### GLU.this.Shader#uniform2F(name, f1, f2)
Provide a vec2 uniform (via uniform2f) to the currently bound shader

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f1 | <code>number</code> | The first float value |
| f2 | <code>number</code> | The second float value |

<a name="GLU.this.Shader+uniform1Fv"></a>

### GLU.this.Shader#uniform1Fv(name, fvec)
Provide an array of floats (via uniform1Fv) to the currently bound shader
  i.e. the shader declares e.g. `uniform float values[19];`

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec | <code>Float32Array</code> | An array of floats |

<a name="GLU.this.Shader+uniform2Fv"></a>

### GLU.this.Shader#uniform2Fv(name, fvec2)
Provide an array of vec2 (via uniform2fv) to the currently bound shader
  i.e. the shader declares e.g. `uniform vec2 vectors[19];`

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec2 | <code>Float32Array</code> | An array of floats, 2 per vector |

<a name="GLU.this.Shader+uniform3F"></a>

### GLU.this.Shader#uniform3F(name, f1, f2, f3)
Provide a vec3 uniform (via uniform3f) to the currently bound shader

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f1 | <code>number</code> | The first float value |
| f2 | <code>number</code> | The second float value |
| f3 | <code>number</code> | The third float value |

<a name="GLU.this.Shader+uniform3Fv"></a>

### GLU.this.Shader#uniform3Fv(name, fvec3)
Provide an array of vec3 (via uniform3fv) to the currently bound shader
  i.e. the shader declares e.g. `uniform vec3 vectors[19];`

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec3 | <code>Float32Array</code> | An array of floats, 3 per vector |

<a name="GLU.this.Shader+uniform4F"></a>

### GLU.this.Shader#uniform4F(name, f1, f2, f3, f4)
Provide a vec4 uniform (via uniform4F) to the currently bound shader

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| f1 | <code>number</code> | The first float value |
| f2 | <code>number</code> | The second float value |
| f3 | <code>number</code> | The third float value |
| f4 | <code>number</code> | The fourth float value |

<a name="GLU.this.Shader+uniform4Fv"></a>

### GLU.this.Shader#uniform4Fv(name, fvec4)
Provide an array of vec4 (via uniform4fv) to the currently bound shader
  i.e. the shader declares e.g. `uniform vec4 vectors[19];`

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| fvec4 | <code>Float32Array</code> | An array of floats, 4 per vector |

<a name="GLU.this.Shader+uniformMatrix4fv"></a>

### GLU.this.Shader#uniformMatrix4fv(name, matrixArray16)
Provide a matrix (via uniformMatrix4fv) to the currently bound shader
 i.e. the shader declares e.g. `uniform mat4 matrix;`

**Kind**: static method of [<code>GLU</code>](#GLU)  

| Param | Type | Description |
| --- | --- | --- |
| name | <code>string</code> | The name of the uniform variable |
| matrixArray16 | <code>Float32Array</code> | An array of 16 floats |


