

function SphereCrystalScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.sphereSpacing = 2.0;
	this._settings.sphereRadius = 1.0;
	this._settings.bulkRadius = 100.0;
	
}

// NB, every function is mandatory and must be defined.

SphereCrystalScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
SphereCrystalScene.prototype.sdf = function()
{
	return `
				uniform float _sphereSpacing;
				uniform float _sphereRadius; 
				uniform float _bulkRadius;            

				float sdSphere(vec3 X, float r)                     
				{                                     
					return length(X) - r;       
				}   

				float sphereLattice(vec3 X, float c)
				{
					vec3 q = mod(X,c) - 0.5*vec3(c,c,c);
					return sdSphere(q, _sphereRadius);
				}   

				float SDF(vec3 X)                     
				{                       
					return opI( sdSphere(X, _bulkRadius), 
								sphereLattice(X, _sphereSpacing) );
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
SphereCrystalScene.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_sphereSpacing", this._settings.sphereSpacing);
	traceProgram.uniformF("_sphereRadius", this._settings.sphereRadius);
	traceProgram.uniformF("_bulkRadius", this._settings.bulkRadius);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
SphereCrystalScene.prototype.getScale = function()
{
	return this._settings.bulkRadius;
}

/*
SphereCrystalScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial cam position default for this scene
SphereCrystalScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-150.300, 0.00000, 0.00000));
	laser.setDirection(new THREE.Vector3(1.00000, 1.66533e-16, -8.67362e-19));
	laser.setEmissionRadius(1.00016);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(-18.0244, -24.9632, -4.79777);
	camera.position.set(-130.310, 66.8365, 81.0769);
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
SphereCrystalScene.prototype.initGui = function(parentFolder)
{
	this.nItem = parentFolder.add(this._settings, 'sphereSpacing', 0.01, 100.0);
	this.rItem = parentFolder.add(this._settings, 'sphereRadius', 0.01, 100.0);
	this.cItem = parentFolder.add(this._settings, 'bulkRadius', 10.0, 100.0);
	this.nItem.onChange( function(value) { snelly.reset(); } );
	this.rItem.onChange( function(value) { snelly.reset(); } );
	this.cItem.onChange( function(value) { snelly.reset(); } );
}

SphereCrystalScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.nItem);
	parentFolder.remove(this.rItem);
	parentFolder.remove(this.cItem);
}









