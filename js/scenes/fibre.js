
function FibreScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.fibreRadius = 2.7;
	this._settings.coilRadius = 1.3;
	this._settings.length = 100.0;
	this._settings.twist = 11.0;
}

FibreScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior is 
// defined by the points with SDF<0.0, with a constant refractive index.
FibreScene.prototype.sdf = function()
{
	return `
				uniform float _fibreRadius;  
				uniform float _coilRadius;            
				uniform float _length; 
				uniform float _twist;                           

				float sdCylinder(vec3 X)                      
				{                                  
					vec2 h = vec2(_fibreRadius, _length);                         
					vec2 d = abs(vec2(length(X.xy), X.z)) - h;         
					return min(max(d.x,d.y),0.0) + length(max(d,0.0)); 
				}       

				float opCoil( vec3 p )
				{
				    float c = cos(_twist*p.z/_length);
				    float s = sin(_twist*p.z/_length);
				    vec3 q = p + 4.0*_coilRadius*vec3(c, s, 0.0);
				    return sdCylinder(q);
				}

				float SDF(vec3 X)                     
				{                           
					return opCoil(X);
				}                                                
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
FibreScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_fibreRadius", this._settings.fibreRadius);
	traceProgram.uniformF("_coilRadius", this._settings.coilRadius);
	traceProgram.uniformF("_length", this._settings.length);
	traceProgram.uniformF("_twist", this._settings.twist);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
FibreScene.prototype.getScale = function()
{
	return Math.min(this._settings.fibreRadius, this._settings.coilRadius);
}

/*
FibreScene.prototype.getBox = function()
{
	var L = this._settings.length;
	var r = this._settings.radius;
	var min = new THREE.Vector3(-3.0*r, -3.0*r, -0.666*L);
	var max = new THREE.Vector3( 3.0*r,  3.0*r,  0.666*L);
	return new THREE.Box3(min, max);
}
*/

// Initial cam position default for this scene
FibreScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(0.272821, -0.782129, -115.617));
	laser.setTarget(new THREE.Vector3(2.27983, -5.68264, -100.000));
	laser.setEmissionRadius(0.69133);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(-23.9941, 9.32417, -130.227);
	camera.position.set(-28.3761, 14.5968, -147.678);
}


// set up gui and callbacks for this scene
FibreScene.prototype.initGui = function(parentFolder)
{
	this.fRadiusItem = parentFolder.add(this._settings, 'fibreRadius', 0.01, 8.0);
	this.fRadiusItem.onChange( function(value) { snelly.reset(); } );

	this.cRadiusItem = parentFolder.add(this._settings, 'coilRadius', 0.01, 2.0);
	this.cRadiusItem.onChange( function(value) { snelly.reset(); } );

	this.lengthItem = parentFolder.add(this._settings, 'length', 0.01, 200.0);
	this.lengthItem.onChange( function(value) { snelly.reset(); } );

	this.twistItem = parentFolder.add(this._settings, 'twist', 0.0, 20.0);
	this.twistItem.onChange( function(value) { snelly.reset(); } );
}

FibreScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.fRadiusItem);
	parentFolder.remove(this.cRadiusItem);
	parentFolder.remove(this.lengthItem);
	parentFolder.remove(this.twistItem);
}








