
We define the rendered scene by specifying, via the {@link Scene#shader}, three GLSL functions:

	- `SDF_SURFACE(vec3 X)`: the SDF of the uber-surface material
	- `SDF_METAL(vec3 X)`: the SDF of the (selected) metal material
	- `SDF_DIELECTRIC(vec3 X)`: the SDF of the (selected) dielectric material
  
These functions are assumed to be SDFs where the negative region corresponds to the interior of the body.
Thus there are at most only three types of material in the scene.

The details of the properties of the three material types can be selected via the {@link Materials} object.
In addition, spatial dependence of the material surface properties can be introduced by providing modulating GLSL functions.

