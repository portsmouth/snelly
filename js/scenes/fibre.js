
function FibreScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.radius = 1.36;
	this._settings.length = 100.0;
	this._settings.twist = 0.0;
}

FibreScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior is 
// defined by the points with SDF<0.0, with a constant refractive index.
FibreScene.prototype.sdf = function()
{
	return `
				uniform float _radius;                 
				uniform float _length;                 

				float SDF(vec3 X)                      
				{                                  
					vec2 h = vec2(_radius, _length);                         
					vec2 d = abs(vec2(length(X.xy), X.z)) - h;         
					return min(max(d.x,d.y),0.0) + length(max(d,0.0)); 
				}                                                      
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
FibreScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_radius", this._settings.radius);
	traceProgram.uniformF("_length", this._settings.length);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
FibreScene.prototype.getScale = function()
{
	return Math.max(this._settings.radius, this._settings.length);
}


FibreScene.prototype.getBox = function()
{
	var L = this._settings.length;
	var r = this._settings.radius;
	var min = new THREE.Vector3(-3.0*r, -3.0*r, -0.666*L);
	var max = new THREE.Vector3( 3.0*r,  3.0*r,  0.666*L);
	return new THREE.Box3(min, max);
}

// Initial cam position default for this scene
FibreScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-6.48839, 4.14743, -114.357));
	laser.setTarget(new THREE.Vector3(0.188796, 0.149697, -100.000));
	laser.setEmissionRadius(1.00000);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(8.08885, -34.1258, -8.56345);
	camera.position.set(-15.6904, 19.2269, -132.803);
}


// set up gui and callbacks for this scene
FibreScene.prototype.initGui = function(parentFolder)
{
	this.radiusItem = parentFolder.add(this._settings, 'radius', 0.01, 10.0);
	this.radiusItem.onChange( function(value) { snelly.reset(); } );

	this.lengthItem = parentFolder.add(this._settings, 'length', 0.01, 100.0);
	this.lengthItem.onChange( function(value) { snelly.reset(); } );
}

FibreScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.radiusItem);
	parentFolder.remove(this.lengthItem);
}







