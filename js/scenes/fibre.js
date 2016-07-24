
function FibreScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.radius = 1.0;
	this._settings.length = 100.0;
}

FibreScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior is 
// defined by the points with SDF<0.0, with a constant refractive index.
FibreScene.prototype.sdf = function()
{
	return `
					uniform float _radius;                 
					uniform float _length;                 
					float SDF(vec3 X)                      
					{                                      
						vec2 h(_radius, _length);                          
						vec2 d = abs(vec2(length(p.xy), p.x)) - h;         
						return min(max(d.x,d.y),0.0) + length(max(d,0.0)); 
					}                                                      
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
FibreScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_radius", this._settings.radius);
	traceProgram.uniformF("_length", this._settings.length);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
FibreScene.prototype.getScale = function()
{
	return this._settings.radius;
}


// Initial cam position default for this scene
FibreScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(50.0, 0.0, 0.0)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
FibreScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(-100.0, 0.0, 0.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
}



// set up gui and callbacks for this scene
FibreScene.prototype.initGui = function(parentFolder)
{

}

// set up gui and callbacks for this material
LinearMetal.prototype.initGui = function(parentFolder)
{
	parentFolder.add(this, 'n400', 0.0, 10.0);
	parentFolder.add(this, 'n700', 0.0, 10.0);
	parentFolder.add(this, 'k400', 0.0, 1000.0);
	parentFolder.add(this, 'k700', 0.0, 1000.0);
}








