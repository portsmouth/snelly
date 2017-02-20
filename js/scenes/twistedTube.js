
function TwistedTubeScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.bounds = new THREE.Vector3(2.5, 2.5, 10.0);
	this._settings.twist = 1.0;
}

// NB, every function is mandatory and must be defined.

TwistedTubeScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
TwistedTubeScene.prototype.sdf = function()
{
	return `
				uniform vec3 _bounds;   
				uniform float _twist;    

				float sdBox(vec3 X, vec3 bounds)                     
				{                                     
					vec3 d = abs(X) - bounds;
					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
				} 

				float opTwist( vec3 p )
				{
				    float c = cos(_twist*p.y/_bounds.y);
				    float s = sin(_twist*p.y/_bounds.y);
				    mat2  m = mat2(c,-s,s,c);
				    vec3  q = vec3(m*p.xz,p.y);
				    return sdBox(q, _bounds);
				}

				float SDF_DIELE(vec3 X)                     
				{                           
					return opTwist(X);
				}      

				float SDF_METAL(vec3 X) { return HUGE_VAL; }
				float SDF_DIFFU(vec3 X) { return HUGE_VAL; }                                  
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
TwistedTubeScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniform3Fv("_bounds", [this._settings.bounds.x, 
										this._settings.bounds.y, 
										this._settings.bounds.z]);
	traceProgram.uniformF("_twist", this._settings.twist);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
TwistedTubeScene.prototype.getScale = function()
{
	var b = this._settings.bounds;
	return Math.max(b.x, b.y, b.z);
}


// Initial cam position default for this scene
TwistedTubeScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(0.109110, 14.3273, 2.11155));
	laser.setTarget(new THREE.Vector3(0.884006, 10.0001, -1.21589));
	laser.setEmissionRadius(0.676533);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(3.30453, 1.84737, 1.72876);
	camera.position.set(-31.8947, 18.7081, 10.8732);
}

/*
TwistedTubeScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// set up gui and callbacks for this scene
TwistedTubeScene.prototype.initGui = function(parentFolder)
{
	this.aItem = parentFolder.add(this._settings.bounds, 'x', 0.01, 10.0);
	this.aItem.onChange( function(value) { snelly.reset(); } );

	this.bItem = parentFolder.add(this._settings.bounds, 'y', 0.01, 10.0);
	this.bItem.onChange( function(value) { snelly.reset(); } );

	this.cItem = parentFolder.add(this._settings.bounds, 'z', 0.01, 10.0);
	this.cItem.onChange( function(value) { snelly.reset(); } );

	this.twistItem = parentFolder.add(this._settings, 'twist', 0.0, 1.0);
	this.twistItem.onChange( function(value) { snelly.reset(); } );
}

TwistedTubeScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.aItem);
	parentFolder.remove(this.bItem);
	parentFolder.remove(this.cItem);
	parentFolder.remove(this.twistItem);
}









