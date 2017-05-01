

function SuperEllipsoidScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.a = 2.5;
	this._settings.b = 1.9;
	this._settings.c = 3.5;
	this._settings.power = 4.5;
}

// NB, every function is mandatory and must be defined.

SuperEllipsoidScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
SuperEllipsoidScene.prototype.sdf = function()
{
	return `
				uniform float _a;
				uniform float _b; 
				uniform float _c;
				uniform float _p;            

				float sdSuperEllipsoid(vec3 X, float a, float b, float c, float p)                     
				{                         
					float xp = pow(abs(X.x/a),p);
					float yp = pow(abs(X.y/b),p);
					float zp = pow(abs(X.z/c),p);     
					return pow(xp + yp + zp, 1.0/p)  - 1.0;       
				}      

				float SDF_DIELE(vec3 X)                     
				{                       
					return sdSuperEllipsoid(X, _a, _b, _c, _p);       
				} 
				
				float SDF_METAL(vec3 X) { return HUGE_VAL; }
				float SDF_DIFFU(vec3 X) { return HUGE_VAL; } //sdBox(X, vec3(-100.0, -10.0, -100.0), vec3(100.0, -1.0, 100.0)); }                                    
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
SuperEllipsoidScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_a", this._settings.a);
	traceProgram.uniformF("_b", this._settings.b);
	traceProgram.uniformF("_c", this._settings.c);
	traceProgram.uniformF("_p", this._settings.power);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
SuperEllipsoidScene.prototype.getScale = function()
{
	return Math.max(this._settings.a, this._settings.b, this._settings.c);
}

/*
SuperEllipsoidScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial cam position default for this scene
SuperEllipsoidScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-8.88858, 0.101274, -0.00490936));
	laser.setTarget(new THREE.Vector3(-2.22387, 0.751606, 2.90733));
	laser.setEmissionRadius(0.00000);
	laser.setEmissionSpreadAngle(0.200000);
	controls.target.set(0.130605, 0.0139136, 0.856111);
	camera.position.set(4.42667, 5.33090, 9.62628);
}


// set up gui and callbacks for this scene
SuperEllipsoidScene.prototype.initGui = function(parentFolder)
{
	this.aItem = parentFolder.add(this._settings, 'a', 0.01, 10.0);
	this.bItem = parentFolder.add(this._settings, 'b', 0.01, 10.0);
	this.cItem = parentFolder.add(this._settings, 'c', 0.01, 10.0);
	this.pItem = parentFolder.add(this._settings, 'power', 2.0, 30.0);
	this.aItem.onChange( function(value) { snelly.reset(); } );
	this.bItem.onChange( function(value) { snelly.reset(); } );
	this.cItem.onChange( function(value) { snelly.reset(); } );
	this.pItem.onChange( function(value) { snelly.reset(); } );
}

SuperEllipsoidScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.aItem);
	parentFolder.remove(this.bItem);
	parentFolder.remove(this.cItem);
	parentFolder.remove(this.pItem);
}









