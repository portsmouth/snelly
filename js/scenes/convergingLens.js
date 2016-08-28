

function ConvergingLensScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.radiusA   = 10.0;
	this._settings.radiusB   = 10.0;
	this._settings.thickness = 1.0;
}

// NB, every function is mandatory and must be defined.

ConvergingLensScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
ConvergingLensScene.prototype.sdf = function()
{
	return `
				uniform float _radiusA;
				uniform float _radiusB;  
				uniform float _thickness;    

 				float sdSphere(vec3 X, float r)
				{
					return length(X) - r;       
				}    

				float SDF(vec3 X)
				{ 
					float t = min(_thickness, _radiusA + _radiusB);
					float l = _radiusA + _radiusB - t;
					float sA = sdSphere(X + vec3( _radiusA - 0.5*t, 0.0, 0.0), _radiusA);
					float sB = sdSphere(X + vec3(-_radiusB + 0.5*t, 0.0, 0.0), _radiusB);
					return opI(sA, sB);
				}
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
ConvergingLensScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_radiusA", this._settings.radiusA);
	traceProgram.uniformF("_radiusB", this._settings.radiusB);
	traceProgram.uniformF("_thickness", this._settings.thickness);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
ConvergingLensScene.prototype.getScale = function()
{
	return Math.max(this._settings.radiusA, this._settings.radiusB);
}


// Initial cam position default for this scene
ConvergingLensScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-13.3410, 0.0, 0.0));
	laser.setDirection(new THREE.Vector3(1.00000, 3.33067e-16, 0.00000));
	laser.setEmissionRadius(3.00000);
	laser.setEmissionSpreadAngle(1.01480);
	controls.target.set(6.23887, -1.60404, 1.39647);
	camera.position.set(18.1117, 9.12804, 33.8416);
}

/*
ConvergingLensScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/



// set up gui and callbacks for this scene
ConvergingLensScene.prototype.initGui = function(parentFolder)
{
	this.radiusAItem   = parentFolder.add(this._settings, 'radiusA', 1.0, 50.0);
	this.radiusBItem   = parentFolder.add(this._settings, 'radiusB', 1.0, 50.0);
	this.thicknessItem = parentFolder.add(this._settings, 'thickness', 0.1, 10.0);

	this.radiusAItem.onChange( function(value) { snelly.reset(); } );
	this.radiusBItem.onChange( function(value) { snelly.reset(); } );
	this.thicknessItem.onChange( function(value) { snelly.reset(); } );
}

ConvergingLensScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.radiusAItem);
	parentFolder.remove(this.radiusBItem);
	parentFolder.remove(this.thicknessItem);
}









