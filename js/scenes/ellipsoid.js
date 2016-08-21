

function EllipsoidScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.a = 1.0;
	this._settings.b = 2.0;
	this._settings.c = 3.0;
}

// NB, every function is mandatory and must be defined.

EllipsoidScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
EllipsoidScene.prototype.sdf = function()
{
	return `
				uniform float _a;
				uniform float _b; 
				uniform float _c;            

				float sphere(vec3 X, float radius)                     
				{                                     
					return length(X) - radius;       
				}      

				float SDF(vec3 X)                     
				{                       
					float x = X.x/_a; float y = X.y/_b; float z = X.z/_c;
					return sphere(vec3(x, y, z), 1.0);       
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
EllipsoidScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_a", this._settings.a);
	traceProgram.uniformF("_b", this._settings.b);
	traceProgram.uniformF("_c", this._settings.c);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
EllipsoidScene.prototype.getScale = function()
{
	return Math.max(this._settings.a, this._settings.b, this._settings.c);
}

/*
EllipsoidScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial cam position default for this scene
EllipsoidScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(-10.0, 10.0, 10.0)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
EllipsoidScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(-6.0, 0.0, 0.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
EllipsoidScene.prototype.initGui = function(parentFolder)
{
	this.aItem = parentFolder.add(this._settings, 'a', 0.01, 10.0);
	this.bItem = parentFolder.add(this._settings, 'b', 0.01, 10.0);
	this.cItem = parentFolder.add(this._settings, 'c', 0.01, 10.0);
	this.aItem.onChange( function(value) { snelly.reset(); } );
	this.bItem.onChange( function(value) { snelly.reset(); } );
	this.cItem.onChange( function(value) { snelly.reset(); } );
}

EllipsoidScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.aItem);
	parentFolder.remove(this.bItem);
	parentFolder.remove(this.cItem);
}









