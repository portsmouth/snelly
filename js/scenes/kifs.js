

function KIFSScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.time = 38.3;
	this._settings.rotAngle = 40.0;
	this._settings.amplitude = 0.6;
	this._settings.iterations = 20;
}

// NB, every function is mandatory and must be defined.

KIFSScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
KIFSScene.prototype.sdf = function()
{
	return `
				// Borrowed from https://www.shadertoy.com/view/Mtf3Rr by Kali 
				// Kaleidoscopic Iterated Function Systems
				uniform float _time;   
				uniform float _rotAngle; 
				uniform float _amplitude; 

				const float Scale=1.25;
				const vec3 Julia=vec3(-3.,-1.5,-0.5);
				const vec3 RotVector=vec3(0.5,-0.05,-1.);
				const float Speed=1.3;

				mat3 rotationMatrix3(vec3 v, float angle)
				{
					float c = cos(radians(angle));
					float s = sin(radians(angle));
					return mat3(c + (1.0 - c)*v.x*v.x,             (1.0 - c)*v.x*v.y - s*v.z,     (1.0 - c)*v.x*v.z + s*v.y,
						            (1.0 - c)*v.y*v.x + s*v.z, c + (1.0 - c)*v.y*v.y,             (1.0 - c)*v.y*v.z - s*v.x,
						            (1.0 - c)*v.z*v.x - s*v.y,     (1.0 - c)*v.z*v.y + s*v.x, c + (1.0 - c)*v.z*v.z
						);
				}

				float SDF_DIELE(vec3 p) 
				{
					p=p.zxy;
					float a=1.5+sin(_time*.3578)*.5;
					p.xy=p.xy*mat2(cos(a),sin(a),-sin(a),cos(a));
					p.x*=.75;
					float time=_time*Speed;
					vec3 ani;
					ani=vec3(sin(time),sin(time),cos(time))*_amplitude;
					p+=sin(p*3.+time*6.)*.04;
					mat3 rot = rotationMatrix3(normalize(RotVector+ani), _rotAngle+sin(time)*10.);
					vec3 pp=p;
					float l;
					const int iter = ${Math.floor(this._settings.iterations)};
					for (int i=0; i<iter; i++) {
						p.xy=abs(p.xy);
						p=p*Scale+Julia;
						p*=rot;
						l=length(p);
					}
					return l*pow(Scale, -float(iter))-.1;
				}            

				float SDF_METAL(vec3 X) { return HUGE_VAL; }
				float SDF_DIFFU(vec3 X) { return sdBox(X, vec3(-100.0, -10.0, -100.0), vec3(100.0, -9.0, 100.0)); }
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
KIFSScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_time", this._settings.time);
	traceProgram.uniformF("_rotAngle", this._settings.rotAngle);
	traceProgram.uniformF("_amplitude", this._settings.amplitude);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
KIFSScene.prototype.getScale = function()
{
	return 10.0;
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
KIFSScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(0.968744, 2.09333, 6.48014));
	laser.setTarget(new THREE.Vector3(-0.313957, -0.694221, 4.04162));
	laser.setEmissionRadius(0.0100000);
	laser.setEmissionSpreadAngle(0.471590);
	controls.target.set(-1.87667, -3.77846, -2.60251);
	camera.position.set(9.93377, 5.37007, 6.78323);
}


// set up gui and callbacks for this scene
KIFSScene.prototype.initGui = function(parentFolder)
{
	this.tItem = parentFolder.add(this._settings, 'time', 0.0, 100.0, 0.01);
	this.tItem.onChange( function(value) { snelly.reset(); } );

	this.rItem = parentFolder.add(this._settings, 'rotAngle', 0.0, 100.0, 0.1);
	this.rItem.onChange( function(value) { snelly.reset(); } );

	this.aItem = parentFolder.add(this._settings, 'amplitude', 0.0, 1.0, 0.01);
	this.aItem.onChange( function(value) { snelly.reset(); } );

	this.iItem = parentFolder.add(this._settings, 'iterations', 1, 100, 1);
	this.iItem.onChange( function(value) { snelly.reset(); } );
}

KIFSScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.tItem);
	parentFolder.remove(this.rItem);
	parentFolder.remove(this.aItem);
	parentFolder.remove(this.iItem);
}









