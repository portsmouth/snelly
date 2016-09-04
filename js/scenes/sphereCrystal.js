

function SphereCrystalScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.sphereSpacing = 2.0;
	this._settings.sphereRadius = 1.0;
	this._settings.width = 9.0;
	this._settings.height = 3.0;
	this._settings.depth = 3.0;
	this._settings.offset = 0.0;
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
				uniform float _width; 
				uniform float _height;
				uniform float _depth; 
				uniform float _offset;          

				float sdBox(vec3 X, vec3 bounds)                     
				{                                     
					vec3 d = abs(X) - bounds;
					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
				} 

				float sdSphere(vec3 X, float r)                     
				{                                     
					return length(X) - r;       
				}   

				float sphereLattice(vec3 X, float c)
				{
					vec3 r = 0.5*vec3(c,c,c);
					vec3 q = mod(X-(1.0-_offset)*r, c) - r;
					return sdSphere(q, _sphereRadius);
				}   

				float SDF(vec3 X)                     
				{
					return opI( sdBox(X, vec3(_width, _height, _depth)), 
								sphereLattice(X, _sphereSpacing) );
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
SphereCrystalScene.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_sphereSpacing", this._settings.sphereSpacing);
	traceProgram.uniformF("_sphereRadius",  this._settings.sphereRadius);
	traceProgram.uniformF("_width",  this._settings.width);
	traceProgram.uniformF("_height", this._settings.height);
	traceProgram.uniformF("_depth",  this._settings.depth);
	traceProgram.uniformF("_offset",  this._settings.offset);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
SphereCrystalScene.prototype.getScale = function()
{
	return Math.max(this._settings.width, this._settings.height, this._settings.depth);
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
	laser.setPosition(new THREE.Vector3(-18.0135, 0.00000, 0.211178));
	laser.setDirection(new THREE.Vector3(1.00000, 1.66533e-16, 4.33681e-19));
	laser.setEmissionRadius(2.00000);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(2.85618, -3.76422, 1.09919);
	camera.position.set(-25.3435, 12.2641, 34.7719);
}


// set up gui and callbacks for this scene
SphereCrystalScene.prototype.initGui = function(parentFolder)
{
	this.nItem = parentFolder.add(this._settings, 'sphereSpacing', 0.1, 100.0);
	this.rItem = parentFolder.add(this._settings, 'sphereRadius', 0.1, 10.0);
	this.wItem = parentFolder.add(this._settings, 'width', 1.0, 100.0);
	this.hItem = parentFolder.add(this._settings, 'height', 1.0, 100.0);
	this.dItem = parentFolder.add(this._settings, 'depth', 1.0, 100.0);
	this.oItem = parentFolder.add(this._settings, 'offset', 0.0, 1.0);
	this.nItem.onChange( function(value) { snelly.reset(); } );
	this.rItem.onChange( function(value) { snelly.reset(); } );
	this.wItem.onChange( function(value) { snelly.reset(); } );
	this.hItem.onChange( function(value) { snelly.reset(); } );
	this.dItem.onChange( function(value) { snelly.reset(); } );
	this.oItem.onChange( function(value) { snelly.reset(); } );
}

SphereCrystalScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.nItem);
	parentFolder.remove(this.rItem);
	parentFolder.remove(this.wItem);
	parentFolder.remove(this.hItem);
	parentFolder.remove(this.dItem);
	parentFolder.remove(this.oItem);
}









