

function WineGlassScene(name, desc) 
{
	Scene.call(this, name, desc);
	this._settings.scaleHeight = 1.0;
}

// NB, every function is mandatory and must be defined.

WineGlassScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
WineGlassScene.prototype.sdf = function()
{
	return `
			uniform float _height;

			float smin (in float a, in float b, in float k) 
			{
				float h = clamp (0.5 + 0.5 * (a - b) / k, 0.0, 1.0);
				return mix (a, b, h) - k * h * (1.0 - h);
			}

			float SDF(vec3 p)
			{
				p.y /= _height;
				float d = length (p);
				float dxz = length (p.xz);
				d = max (max (d - 1.45, 1.4 - d), p.y - 0.5);
				d = min (d, max (dxz - 1.0, p.y + 3.0));
				d = smin (d, max (dxz - 0.2, p.y + 1.4), 0.2);
				return max (d, -p.y - 3.05);
			}                                    
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
WineGlassScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_height", this._settings.scaleHeight);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
WineGlassScene.prototype.getScale = function()
{
	return this._settings.scaleHeight;
}


WineGlassScene.prototype.getBox = function()
{
	var h = this._settings.scaleHeight;
	var min = new THREE.Vector3(-1.5*h, -3.5*h, -1.5*h);
	var max = new THREE.Vector3(1.5*h, 1.5*h, 1.5*h);
	return new THREE.Box3(min, max);
}


// Initial cam position default for this scene
WineGlassScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-1.65450, 1.69287, -1.90594));
	laser.setDirection(new THREE.Vector3(0.400671, -0.768580, 0.498746));
	laser.setEmissionRadius(0.00000);
	laser.setEmissionSpreadAngle(2.00000);
	controls.target.set(1.30918, -1.19405, 0.224590);
	camera.position.set(-8.30344, 2.52862, 2.53246);
}


// set up gui and callbacks for this scene
WineGlassScene.prototype.initGui = function(parentFolder)
{
	this.heightItem = parentFolder.add(this._settings, 'scaleHeight', 1.0, 2.0);
	this.heightItem.onChange( function(value) { snelly.reset(); } );
}

WineGlassScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.heightItem);
}









