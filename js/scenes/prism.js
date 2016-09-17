


function PrismScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.sides = 3;
	this._settings.radius = 0.6;
	this._settings.height = 10;
	this._settings.twist = 4;
}

// NB, every function is mandatory and must be defined.

PrismScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
PrismScene.prototype.sdf = function()
{
	return `
				uniform float _radius;
				uniform float _height;
				uniform float _twist;

				vec3 rotate2d(in vec3 X, in float angle) 
				{
					float c = cos(angle);
					float s = sin(angle);
					return vec3(c*X.x - s*X.y, c*X.y + s*X.x, X.z);
				}

				float sdWedge(vec3 X, float phi)
				{
					const float eps = 1.0e-4;
					float sh = max(X.z-_height, -X.z-_height);
					float s1 = -X.x-eps;
					float s2 = X.x*cos(phi) - X.y*sin(phi) - eps;
					float s3 = X.x*sin(0.5*phi) + (X.y-_radius)*cos(0.5*phi) - eps;
					return max(max(max(s1, s2), s3), sh);
				}

				float sdPrism(vec3 X)                    
				{               
					const int sides = ${Math.floor(this._settings.sides)};
					float phi = 2.0*M_PI/float(sides);
					float s = sdWedge(X, phi);
					for (int n=1; n<sides; n++) 
					{
						float w = sdWedge(rotate2d(X, float(n)*phi), phi);
						s = opU(s, w);
					}
					return s;
				}      

				float opTwist( vec3 p )
				{
				    float c = cos(_twist*p.z/_height);
				    float s = sin(_twist*p.z/_height);
				    mat2  m = mat2(c,-s,s,c);
    				vec3  q = vec3(m*p.xy,p.z);
				    return sdPrism(q);
				}

				float SDF(vec3 X)                     
				{                           
					return opTwist(X);
				}   

                               
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
PrismScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_radius", this._settings.radius);
	traceProgram.uniformF("_height", this._settings.height);
	traceProgram.uniformF("_twist", this._settings.twist);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
PrismScene.prototype.getScale = function()
{
	return Math.max(this._settings.radius, this._settings.height);
}




// Initial laser position and direction defaults for this scene
PrismScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-22.5525, 15.4152, 23.1324));
	laser.setTarget(new THREE.Vector3(-0.323230, -0.302321, 10.0001));
	laser.setEmissionRadius(0.00000);
	laser.setEmissionSpreadAngle(0.200000);
	controls.target.set(-6.16008, -1.63005, 4.21104);
	camera.position.set(1.20902, 0.422143, -11.6874);
}


// set up gui and callbacks for this scene
PrismScene.prototype.initGui = function(parentFolder)
{
	this.sidesItem = parentFolder.add(this._settings, 'sides', 3, 20);
	this.sidesItem.onChange( function(value) { snelly.reset(); } );

	this.radiusItem = parentFolder.add(this._settings, 'radius', 0.1, 2.0);
	this.radiusItem.onChange( function(value) { snelly.reset(); } );

	this.heightItem = parentFolder.add(this._settings, 'height', 0.0, 10.0);
	this.heightItem.onChange( function(value) { snelly.reset(); } );

	this.twistItem = parentFolder.add(this._settings, 'twist', 0.0, 20.0);
	this.twistItem.onChange( function(value) { snelly.reset(); } );
}

PrismScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.sidesItem);
	parentFolder.remove(this.radiusItem);
	parentFolder.remove(this.heightItem);
	parentFolder.remove(this.twistItem);
}














