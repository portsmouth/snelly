
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


// Initial cam position default for this scene
FibreScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(-13.0, 5.0, -39.0)
	controls.target.set(-4.0, -1.0, -16.0);
}


// Initial laser position and direction defaults for this scene
FibreScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(0.0, -1.0, -24.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));

	Scene.prototype.setLaser.call(this, laser);
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







