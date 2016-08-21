

function StackScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.extent     = 10.0;
	this._settings.thickness  = 1.0;
	this._settings.separation = 1.0;
	this._settings.NUM_LAYERS = 6;
}

// NB, every function is mandatory and must be defined.

StackScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
StackScene.prototype.sdf = function()
{
	return `
				uniform float _extent;   
				uniform float _thickness;
				uniform float _separation;

				float sdBox(vec3 X, vec3 bounds)                     
				{                                     
					vec3 d = abs(X) - bounds;
					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
				}              

				float SDF(vec3 X)                     
				{               
					float inf = 1.0e9;
					float sd = inf;
					vec3 bounds = vec3(_extent, 0.5*_thickness, _extent);
					float shift = _thickness + _separation;
					for(int i=0; i < ${Math.floor(this._settings.NUM_LAYERS)}; i++) 
			    	{
			    		vec3 disp = vec3(0.0, float(i)*shift, 0.0);
			    		float slab = sdBox(X+disp, bounds);
			    		sd = min(sd, slab);
			    	}	    
					return sd;       
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
StackScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_extent",     this._settings.extent);
	traceProgram.uniformF("_thickness",  this._settings.thickness);
	traceProgram.uniformF("_separation", this._settings.separation);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
StackScene.prototype.getScale = function()
{
	var w = this._settings.extent;
	var t = this._settings.thickness;
	var s = this._settings.separation;
	var n = this._settings.NUM_LAYERS;
	return Math.max(w, (t + s)*n);
}

/*
StackScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/



// Initial cam position default for this scene
StackScene.prototype.setCam = function(controls, camera)
{
	var s = this.getScale();
	camera.position.set(-1.5*s, 1.5*s, 1.5*s)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
StackScene.prototype.setLaser = function(laser)
{
	var w = this._settings.extent;
	var t = this._settings.thickness;
	var s = this._settings.separation;
	var n = this._settings.NUM_LAYERS;

	laser.setPosition(new THREE.Vector3(0.0, 0.5*w, 0.0));
	laser.setDirection(new THREE.Vector3(0.0, -1.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
StackScene.prototype.initGui = function(parentFolder)
{
	this.extentItem     = parentFolder.add(this._settings, 'extent', 1.0, 30.0);
	this.thicknessItem  = parentFolder.add(this._settings, 'thickness', 0.0, 1.0);
	this.separationItem = parentFolder.add(this._settings, 'separation', 0.0, 1.0);
	this.NUM_LAYERSItem = parentFolder.add(this._settings, 'NUM_LAYERS', 1, 100, 1);

	this.extentItem.onChange( function(value) { snelly.reset(); } );
	this.thicknessItem.onChange( function(value) { snelly.reset(); } );
	this.separationItem.onChange( function(value) { snelly.reset(); } );
	this.NUM_LAYERSItem.onChange( function(value) { snelly.reset(); } );
}

StackScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.extentItem);
	parentFolder.remove(this.thicknessItem);
	parentFolder.remove(this.separationItem);
	parentFolder.remove(this.NUM_LAYERSItem);
}









