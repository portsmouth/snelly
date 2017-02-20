

function MetaShapesScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.time = 0.0;
	this._settings.blend = 1.85;
}

// NB, every function is mandatory and must be defined.

MetaShapesScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
MetaShapesScene.prototype.sdf = function()
{
	return `
				// Borrowed from https://www.shadertoy.com/view/Xls3R7 by LeWiZ 
				uniform float _time;
				uniform float _blend;

				float sphere(vec3 pos)
				{
					return length(pos)-1.0;   
				}

				float box(vec3 pos)
				{
				    vec3 d = abs(pos) - 1.0;
				  	return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
				}

				float torus(vec3 pos)
				{
					vec2 q = vec2(length(pos.xz)-2.0,pos.y);
				  	return length(q)-0.5;   
				}

				float blob3(float d1, float d2, float d3)
				{
				    float k = _blend;
					return 0.5 - (exp(-k*d1)+exp(-k*d2)+exp(-k*d3))/k;
				}

				float SDF_DIELE(vec3 pos)                    
				{                        
      				float t = _time;
    
				    float p = torus(pos + vec3(0.0,3.0,0.0));
					float b = sphere(0.5*(pos + vec3(cos(t*0.5),sin(t*0.3),0.0)));
				    float s = box(2.0*(pos + 3.0 * vec3(cos(t*1.1),cos(t*1.3),cos(t*1.7))))/2.0;

				    return blob3(p, b, s);
				} 
				
				float SDF_METAL(vec3 X) { return HUGE_VAL; }
				float SDF_DIFFU(vec3 X) { return HUGE_VAL; }
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
MetaShapesScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_time", this._settings.time);
	traceProgram.uniformF("_blend", this._settings.blend);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
MetaShapesScene.prototype.getScale = function()
{
	return 10.0;
}


/*
MetaShapesScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial laser position and direction defaults for this scene
MetaShapesScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-5.44442, -1.57662, -3.23921));
	laser.setTarget(new THREE.Vector3(-1.67982, -2.34401, -2.12558));
	laser.setEmissionRadius(0.00000);
	laser.setEmissionSpreadAngle(0.500000);
	controls.target.set(-1.88619, -1.29066, 0.156162);
	camera.position.set(9.37386, 2.10021, -8.10325);
}


// set up gui and callbacks for this scene
MetaShapesScene.prototype.initGui = function(parentFolder)
{
	this.timeItem = parentFolder.add(this._settings, 'time', 0.0, 10.0);
	this.timeItem.onChange( function(value) { snelly.reset(); } );

	this.blendItem = parentFolder.add(this._settings, 'blend', 0.0, 10.0);
	this.blendItem.onChange( function(value) { snelly.reset(); } );
}

MetaShapesScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.timeItem);
	parentFolder.remove(this.blendItem);
}









