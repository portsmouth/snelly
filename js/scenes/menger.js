
function MengerScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.MENGER_ITERATIONS = 4;
	this._settings.tile = false;
	this._settings.tileScale = 1.0;
}

// NB, every function is mandatory and must be defined.

MengerScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
MengerScene.prototype.sdf = function()
{
	var code = `
				float maxcomp(in vec3 p) { return max(p.x,max(p.y,p.z));}

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
					float scale = 1.0;
					const int iter = ${Math.floor(this._settings.MENGER_ITERATIONS)};
					for (int i=0; i<iter; i++) 
					{
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
MengerScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(-3.0, 3.0, 3.0)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
MengerScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(-2.0, 0.0, 0.0));
	laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
MengerScene.prototype.initGui = function(parentFolder)
{
	this.iterItem = parentFolder.add(this._settings, 'MENGER_ITERATIONS', 1, 8, 1);
	this.iterItem.onChange( function(value) { snelly.reset(); } );

	this.tileItem = parentFolder.add(this._settings, 'tile', false);
	this.tileItem.onChange( function(value) { snelly.reset(); } );

	this.tileScaleItem = parentFolder.add(this._settings, 'tileScale', 0.1, 10.0);
	this.tileScaleItem.onChange( function(value) { snelly.reset(); } );
}

MengerScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.iterItem);
	parentFolder.remove(this.tileItem);
	parentFolder.remove(this.tileScaleItem);
}









