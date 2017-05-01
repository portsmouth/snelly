

function GemScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.angularFacets = 8;
	this._settings.scaleWidth = 1.0;
	this._settings.scaleHeight = 1.0;
}

// NB, every function is mandatory and must be defined.

GemScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
GemScene.prototype.sdf = function()
{
	return `
				// Borrowed from https://www.shadertoy.com/view/ltfXDM by Nrx 
				uniform float _scaleWidth;                 
				uniform float _scaleHeight;  

				vec3 vRotateY (in vec3 p, in float angle) 
				{
					float c = cos(angle);
					float s = sin(angle);
					return vec3 (c * p.x - s * p.z, p.y, c * p.z + s * p.x);
				}

				vec3 normalTopA = normalize (vec3 (0.0, 1.0, 1.4));
				vec3 normalTopB = normalize (vec3 (0.0, 1.0, 1.0));
				vec3 normalTopC = normalize (vec3 (0.0, 1.0, 0.5));
				vec3 normalBottomA = normalize (vec3 (0.0, -1.0, 1.0));
				vec3 normalBottomB = normalize (vec3 (0.0, -1.0, 1.6));

				float SDF_METAL(vec3 p)                    
				{            
				    p.xz /= _scaleWidth;
				    p.y  /= _scaleHeight;
					float topCut = p.y - 1.0;
					float angleStep = M_PI / float(${this._settings.angularFacets});
					float angle = angleStep * (0.5 + floor (atan (p.x, p.z) / angleStep));
					vec3 q = vRotateY (p, angle);
					float topA = dot (q, normalTopA) - 2.0;
					float topC = dot (q, normalTopC) - 1.5;
					float bottomA = dot (q, normalBottomA) - 1.7;
					q = vRotateY (p, -angleStep * 0.5);
					angle = angleStep * floor (atan (q.x, q.z) / angleStep);
					q = vRotateY (p, angle);
					float topB = dot (q, normalTopB) - 1.85;
					float bottomB = dot (q, normalBottomB) - 1.9;

					//float box = sdBox(p, vec3(-100.0, -2.5, -100.0), vec3(100.0, -2.0, 100.0));
					float gem  = max(topCut, max(topA, max(topB, max(topC, max (bottomA, bottomB)))));
					return gem;
					//return opU(box, gem);
				}     
				
				float SDF_DIELE(vec3 X) { return HUGE_VAL; }
				float SDF_DIFFU(vec3 X) { return HUGE_VAL; } //sdBox(X, vec3(-100.0, -2.5, -100.0), vec3(100.0, -2.0, 100.0)); }
                         
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
GemScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_scaleWidth", this._settings.scaleWidth);
	traceProgram.uniformF("_scaleHeight", this._settings.scaleHeight);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
GemScene.prototype.getScale = function()
{
	return Math.min(this._settings.scaleHeight, this._settings.scaleWidth);
}


/*
GemScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial laser position and direction defaults for this scene
GemScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(2.14233e-7, 2.27653, -1.42807));
	laser.setTarget(new THREE.Vector3(0.110644, 0.254047, 2.29894));
	laser.setEmissionRadius(0.0100000);
	laser.setEmissionSpreadAngle(1.00000);
	controls.target.set(0.661976, -0.511723, 0.693621);
	camera.position.set(-7.73283, 3.69481, 7.68164);
}


// set up gui and callbacks for this scene
GemScene.prototype.initGui = function(parentFolder)
{
	this.widthItem = parentFolder.add(this._settings, 'scaleWidth', 1.0, 2.0);
	this.heightItem = parentFolder.add(this._settings, 'scaleHeight', 1.0, 2.0);
	this.angularItem = parentFolder.add(this._settings, 'angularFacets', 1, 100, 1);

	this.widthItem.onChange( function(value) { snelly.reset(); } );
	this.heightItem.onChange( function(value) { snelly.reset(); } );
	this.angularItem.onChange( function(value) { snelly.reset(); } );
}

GemScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.widthItem);
	parentFolder.remove(this.heightItem);
	parentFolder.remove(this.angularItem);
}









