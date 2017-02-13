

function TumblerScene(name, desc) 
{
	Scene.call(this, name, desc);
	this._settings.scaleHeight = 1.0;
}

// NB, every function is mandatory and must be defined.

TumblerScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
TumblerScene.prototype.sdf = function()
{
	return `
			// Borrowed from https://www.shadertoy.com/view/4s2GDV by mu6k
			uniform float _height;
			
			float SDF_DIELE(vec3 p)
			{
				p.y /= _height;
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

			float SDF_METAL(vec3 X) { return 1.0; }
			float SDF_DIFFU(vec3 X) { return 1.0; }                       
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
TumblerScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_height", this._settings.scaleHeight);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
TumblerScene.prototype.getScale = function()
{
	return this._settings.scaleHeight;
}

/*
TumblerScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/

// Initial cam position default for this scene
TumblerScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-5.48136, 3.41904, -0.0665851));
	laser.setTarget(new THREE.Vector3(-0.598250, 1.00000, -0.956469));
	laser.setEmissionRadius(0.0100000);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(2.12505, -1.61120, 1.04772);
	camera.position.set(-8.34300, 6.04460, 1.13916);
}


// set up gui and callbacks for this scene
TumblerScene.prototype.initGui = function(parentFolder)
{
	this.heightItem   = parentFolder.add(this._settings, 'scaleHeight', 1.0, 2.0);
	this.heightItem.onChange( function(value) { snelly.reset(); } );
}

TumblerScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.heightItem);
}









