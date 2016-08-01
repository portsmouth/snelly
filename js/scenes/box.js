
function BoxScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.bounds = new THREE.Vector3(1.0, 1.0, 1.0);
}

// NB, every function is mandatory and must be defined.

BoxScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
BoxScene.prototype.sdf = function()
{
	return `
				uniform vec3 _bounds;                

				float SDF(vec3 X)                     
				{                                     
					vec3 d = abs(X) - _bounds;
  					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
BoxScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniform3Fv("_bounds", [this._settings.bounds.x, 
										this._settings.bounds.y, 
										this._settings.bounds.z]);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
BoxScene.prototype.getScale = function()
{
	var b = this._settings.bounds;
	return Math.max(b.x, b.y, b.z);
}


// Initial cam position default for this scene
BoxScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(-5.0, 5.0, 5.0)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
BoxScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(-6.0, 0.0, 0.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
BoxScene.prototype.initGui = function(parentFolder)
{
	this.aItem = parentFolder.add(this._settings.bounds, 'x', 0.01, 10.0);
	this.aItem.onChange( function(value) { renderer.reset(); } );

	this.bItem = parentFolder.add(this._settings.bounds, 'y', 0.01, 10.0);
	this.bItem.onChange( function(value) { renderer.reset(); } );

	this.cItem = parentFolder.add(this._settings.bounds, 'z', 0.01, 10.0);
	this.cItem.onChange( function(value) { renderer.reset(); } );
}

BoxScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.aItem);
	parentFolder.remove(this.bItem);
	parentFolder.remove(this.cItem);
}









