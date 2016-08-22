

function GlassScene(name, desc) 
{
	Scene.call(this, name, desc);
}

// NB, every function is mandatory and must be defined.

GlassScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
GlassScene.prototype.sdf = function()
{
	return `
				float SDF(vec3 p)
				{
					float a = (length(p.xz)-1.0-p.y*.15)*.85;
					a = max(abs(p.y)-1.0,a);
					float a2 = (length(p.xz)-0.9-p.y*.15)*.85;
					a = max(a,-max(-.8-p.y,a2));
					a = max(a,-length(p+vec3(.0,4.0,.0))+3.09);
					
					vec3 p2 = p; p2.xz*=(1.0-p.y*.15);
					float angle = atan(p2.x,p2.z);
					float mag = length(p2.xz);
					angle = mod(angle,3.14159*.125)-3.14159*.125*.5;
					p2.xz = vec2(cos(angle),sin(angle))*mag;
					a = max(a,(-length(p2+vec3(-7.0,0.0,0.0))+6.05)*.85);
					return a;
				}                                  
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
GlassScene.prototype.syncShader = function(traceProgram)
{

}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
GlassScene.prototype.getScale = function()
{
	return 10.0;
}

/*
GlassScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/

// Initial cam position default for this scene
GlassScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(-10.0, 10.0, 10.0)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
WineGlassScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(-6.0, 6.0, -6.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
GlassScene.prototype.initGui = function(parentFolder)
{

}

GlassScene.prototype.eraseGui = function(parentFolder)
{

}









