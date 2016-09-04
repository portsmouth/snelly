

function EllipsoidScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.a = 2.0;
	this._settings.b = 3.0;
	this._settings.c = 5.0;
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

				float sdSphere(vec3 X, float radius)                     
				{                                     
					return length(X) - radius;       
				}      

				float SDF(vec3 X)                     
				{                       
					return sdSphere(vec3(X.x/_a, X.y/_b, X.z/_c), 1.0);       
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
EllipsoidScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-8.88858, 0.101274, -0.00490936));
	laser.setTarget(new THREE.Vector3(0.252166, -2.76093, 1.85217));
	laser.setEmissionRadius(0.250000);
	laser.setEmissionSpreadAngle(1.01480);
	controls.target.set(1.66261, -1.28846, 0.929629);
	camera.position.set(4.60906, -2.68147, 17.3705);
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









