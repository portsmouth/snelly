

function KnotScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.windingNumber = 2.5;
	this._settings.radius = 0.8;
}

// NB, every function is mandatory and must be defined.

KnotScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
KnotScene.prototype.sdf = function()
{
	return `
				// Borrowed from https://www.shadertoy.com/view/XlSGWR by vgs 
				uniform float _k;   
				uniform float _r;             
				#define TWOPI 6.28318530718

				float SDF(vec3 p)                    
				{                        
				    float r = length(p.xy);
				    float oa, a = atan(p.y, p.x); oa = _k*a;
				    a = mod(a, 0.001*TWOPI) - 0.001*TWOPI/2.0;
				    p.xy = r*vec2(cos(a), sin(a)); p.x -= 6.0;
				    p.xz = cos(oa)*p.xz + sin(oa)*vec2(-p.z, p.x);
				    p.x = abs(p.x) - 1.35; 
				    return length(p) - _r;
				}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
KnotScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_k", this._settings.windingNumber);
	traceProgram.uniformF("_r", this._settings.radius);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
KnotScene.prototype.getScale = function()
{
	return this._settings.radius;
}


/*
KnotScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial laser position and direction defaults for this scene
KnotScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-0.412294, -7.09606, 1.35382));
	laser.setDirection(new THREE.Vector3(-0.807465, 0.484682, -0.336279));
	laser.setEmissionRadius(0.0100000);
	laser.setEmissionSpreadAngle(0.300000);
	controls.target.set(2.90170, -0.573916, -0.616009);
	camera.position.set(4.95500, -0.404701, 23.2964);
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
KnotScene.prototype.initGui = function(parentFolder)
{
	this.kItem = parentFolder.add(this._settings, 'windingNumber', 1.0, 6.0, 0.5);
	this.rItem = parentFolder.add(this._settings, 'radius', 0.1, 6.0);
	this.kItem.onChange( function(value) { snelly.reset(); } );
	this.rItem.onChange( function(value) { snelly.reset(); } );
}

KnotScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.kItem);
	parentFolder.remove(this.rItem);
}









