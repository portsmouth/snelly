

function DivergingLensScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.radius   = 10.0;
	this._settings.thickness = 1.0;
}

// NB, every function is mandatory and must be defined.

DivergingLensScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
DivergingLensScene.prototype.sdf = function()
{
	return `
			uniform float _radius;
			uniform float _thickness;    
			
			float sdBox(vec3 X, vec3 bounds)                     
			{                                     
				vec3 d = abs(X) - bounds;
				return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
			} 

				float sdSphere(vec3 X, float r)                     
			{                                     
				return length(X) - r;       
			}    

			float SDF(vec3 X)                     
			{            
				float r = max(_radius, 3.0*_thickness);
				float y = 0.333*r;
				vec3 b = vec3(0.25*r, y, y);
				float block = sdBox(X, b);
				float sA = sdSphere(X + vec3( r + 0.5*_thickness, 0.0, 0.0), r); 
				float sB = sdSphere(X + vec3(-r - 0.5*_thickness, 0.0, 0.0), r);   
				return opS(block, opU(sA, sB));                       
			}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
DivergingLensScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_radius", this._settings.radius);
	traceProgram.uniformF("_thickness", this._settings.thickness);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
DivergingLensScene.prototype.getScale = function()
{
	return this._settings.radius;
}


// Initial cam position default for this scene
DivergingLensScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-13.3410, 0.101274, -0.00490936));
	laser.setDirection(new THREE.Vector3(1.00000, 3.33067e-16, 0.00000));
	laser.setEmissionRadius(3.00000);
	laser.setEmissionSpreadAngle(1.01480);
	controls.target.set(6.23887, -1.60404, 1.39647);
	camera.position.set(18.1117, 9.12804, 33.8416);
}

/*
DivergingLensScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/

// set up gui and callbacks for this scene
DivergingLensScene.prototype.initGui = function(parentFolder)
{
	this.radiusItem   = parentFolder.add(this._settings, 'radius', 1.0, 50.0);
	this.thicknessItem = parentFolder.add(this._settings, 'thickness', 0.0, 10.0);

	this.radiusItem.onChange( function(value) { snelly.reset(); } );
	this.thicknessItem.onChange( function(value) { snelly.reset(); } );
}

DivergingLensScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.radiusItem);
	parentFolder.remove(this.thicknessItem);
}









