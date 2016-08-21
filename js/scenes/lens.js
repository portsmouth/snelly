

function LensScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.radius = 5.0;
}

// NB, every function is mandatory and must be defined.

LensScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
LensScene.prototype.sdf = function()
{
	return `
				uniform float _radiusA;
				uniform float _radiusB;      

 				float sphere(vec3 X)                     
				{                                     
					return length(X) - _radius;       
				}    

				float SDF(vec3 X)                     
				{                                     
					return length(X) - _radius;       
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
LensScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_radius", this._settings.radius);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
LensScene.prototype.getScale = function()
{
	return this._settings.radius;
}


// Initial cam position default for this scene
LensScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(-10.0, 10.0, 10.0)
	controls.target.set(0.0, 0.0, 0.0);
}

/*
LensScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial laser position and direction defaults for this scene
LensScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(-6.0, 0.0, 0.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
LensScene.prototype.initGui = function(parentFolder)
{
	this.radiusItem = parentFolder.add(this._settings, 'radius', 0.01, 10.0);
	this.radiusItem.onChange( function(value) { snelly.reset(); } );
}

LensScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.radiusItem);
}









