
function MengerScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.iterations = 4;
	this._settings.tile = false;
	this._settings.tileScale = 1.0;
	this._settings.rotate = 0.0;
	this._settings.scale = 1.0;
}

// NB, every function is mandatory and must be defined.

MengerScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
MengerScene.prototype.sdf = function()
{
	var code = `
				uniform float _rotate;
				uniform float _scale;
				#define TWOPI 6.28318530718

				vec2 rotate(vec2 v, float a) {
					return vec2(cos(a)*v.x + sin(a)*v.y, -sin(a)*v.x + cos(a)*v.y); 
				}

				float sdBox(vec3 X, vec3 bounds)                     
				{                                     
					vec3 d = abs(X) - bounds;
					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
				} 

				float sdCross( in vec3 p )
				{
					float inf = 1.0e6;
					float da = sdBox(p.xyz, vec3(inf,1.0,1.0));
					float db = sdBox(p.yzx, vec3(1.0,inf,1.0));
					float dc = sdBox(p.zxy, vec3(1.0,1.0,inf));
					return min(da,min(db,dc));
				}

				float menger(vec3 X) 
				{
					float sd = sdBox(X, vec3(1.0));
					float scale = _scale;
					const int iter = ${Math.floor(this._settings.iterations)};
					for (int i=0; i<iter; i++) 
					{
						X.xy = rotate(X.xy, 2.0*sin(TWOPI*_rotate));

						vec3 a = mod(X*scale, 2.0)-1.0;
      					scale *= 3.0;
      					vec3 r = abs(1.0 - 3.0*abs(a));
						float da = max(r.x,r.y);
						float db = max(r.y,r.z);
						float dc = max(r.z,r.x);
						float c = (min(da,min(db,dc))-1.0)/scale;

						sd = max(sd, c);
					}
					return sd;
				}                                  
	`;

	var sdfCode;
	if (!this._settings.tile)
	{
		sdfCode = code + `
			float SDF(vec3 X)
			{
				return menger(X);
			}
		`;
	}
	else
	{
		sdfCode = code + `
			float opRep( vec3 p, vec3 c )
			{
				vec3 q = mod(p,c)-0.5*c;
				return menger(q);
			}
			float SDF(vec3 X)
			{
				float s = float(${this._settings.tileScale});
				return opRep(X, vec3(s, s, s));
			}
		`;
	}
	return sdfCode;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
MengerScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_rotate", this._settings.rotate);
	traceProgram.uniformF("_scale", this._settings.scale);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
MengerScene.prototype.getScale = function()
{
	return 1.0;
}


/*
MengerScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/



// Initial cam position default for this scene
MengerScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-1.74246, 1.06872, -1.18454));
	laser.setDirection(new THREE.Vector3(0.505430, -0.479445, 0.717407));
	laser.setEmissionRadius(0.0100000);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(2.01584, -0.603169, -0.808580);
	camera.position.set(-4.81888, 2.39746, 3.47067);
}



// set up gui and callbacks for this scene
MengerScene.prototype.initGui = function(parentFolder)
{
	this.iterItem = parentFolder.add(this._settings, 'iterations', 1, 8, 1);
	this.iterItem.onChange( function(value) { snelly.reset(); } );

	this.rotateItem = parentFolder.add(this._settings, 'rotate', 0.0, 1.0);
	this.rotateItem.onChange( function(value) { snelly.reset(); } );

	this.scaleItem = parentFolder.add(this._settings, 'scale', 0.0, 10.0);
	this.scaleItem.onChange( function(value) { snelly.reset(); } );

	this.tileItem = parentFolder.add(this._settings, 'tile', false);
	this.tileItem.onChange( function(value) { snelly.reset(); } );

	this.tileScaleItem = parentFolder.add(this._settings, 'tileScale', 0.1, 10.0);
	this.tileScaleItem.onChange( function(value) { snelly.reset(); } );
}

MengerScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.iterItem);
	parentFolder.remove(this.rotateItem);
	parentFolder.remove(this.scaleItem);
	parentFolder.remove(this.tileItem);
	parentFolder.remove(this.tileScaleItem);
}









