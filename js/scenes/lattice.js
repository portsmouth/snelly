

function LatticeScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.shapes = ['Sphere', 'Cube', 'Octahedron', 'Dodecahedron', 'Icosahedron', 
								'TruncatedOctahedron', 'TruncatedIcosahedron'];
	this._settings.shape = 'Sphere';
	this._settings.size = 1.0;
	this._settings.spacing = 2.0;
	this._settings.width = 3.0;
	this._settings.height = 7.0;
	this._settings.depth = 9.0;
	this._settings.offset = 0.0;
	this._settings.bulge = 0.25;
}

// NB, every function is mandatory and must be defined.

LatticeScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
LatticeScene.prototype.sdf = function()
{
	var func = '';
	var ret = '';
	if (this._settings.bulge < 1.0e-6)
	{
		func = `d = max(d, abs(dot(p, GDFVectors[i])))`;
		ret = `d - r`;
	}
	else
	{
		this._settings.bulge = Math.max(this._settings.bulge, 0.05);
		func = `d += pow(abs(dot(p, GDFVectors[i])), 1.0/e)`;
		ret = `pow(d, e) - r`;

		//func = `d += pow(abs(dot(p, GDFVectors[i])), e)`;
		//ret = `pow(d, 1.0/e) - r`;
	}

	return  `   uniform float _size;
				uniform float _spacing;
				uniform float _width; 
				uniform float _height;
				uniform float _depth; 
				uniform float _offset;
				uniform float _bulge;        

				// Uses code from HG_SDF library, http://mercury.sexy/hg_sdf/ 
				uniform vec3 GDFVectors[19];

				float sdOctahedron(vec3 p, float r, float e) {
					float d = 0.0;
					for (int i = 3; i <= 6; ++i) ${func};
					return ${ret};
				}

				float sdDodecahedron(vec3 p, float r, float e) {
					float d = 0.0;
					for (int i = 13; i <= 18; ++i) ${func};
					return ${ret};
				}

				float sdIcosahedron(vec3 p, float r, float e) {
					float d = 0.0;
					for (int i = 3; i <= 12; ++i) ${func};
					return ${ret};
				}

				float sdTruncatedOctahedron(vec3 p, float r, float e) {
					float d = 0.0;
					for (int i = 0; i <= 6; ++i) ${func};
					return ${ret};
				}

				float sdTruncatedIcosahedron(vec3 p, float r, float e) {
					float d = 0.0;
					for (int i = 3; i <= 18; ++i) ${func};
					return ${ret};
				}

				float sdSphere(vec3 X, float r, float e)                  
				{                                     
					return length(X) - r;       
				}      

				float sdCube(vec3 X, float r, float e)               
				{                                   
					vec3 d = abs(X) - vec3(r,r,r);
					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));           
				}  

				float lattice(vec3 X, float c)
				{
					vec3 r = 0.5*vec3(c,c,c);
					vec3 q = mod(X-(1.0-_offset)*r, c) - r;
					return sd${this._settings.shape}(q, _size, _bulge/5.0);
				}

				float sdBox(vec3 X, vec3 bounds)                     
				{                                     
					vec3 d = abs(X) - bounds;
					return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
				} 

				float SDF(vec3 X)                     
				{
					return opI( sdBox(X, vec3(_width, _height, _depth)), 
								lattice(X, _spacing) );
				}         
			`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
LatticeScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)

	var PHI = Math.sqrt(5.0)*0.5 + 0.5;

	var GDFVectors = [];
	GDFVectors.push(new THREE.Vector3(1, 0, 0).normalize());
	GDFVectors.push(new THREE.Vector3(0, 1, 0).normalize());
	GDFVectors.push(new THREE.Vector3(0, 0, 1).normalize());
	GDFVectors.push(new THREE.Vector3(1, 1, 1).normalize());
	GDFVectors.push(new THREE.Vector3(-1, 1, 1).normalize());
	GDFVectors.push(new THREE.Vector3(1, -1, 1).normalize());
	GDFVectors.push(new THREE.Vector3(1, 1, -1).normalize());
	GDFVectors.push(new THREE.Vector3(0, 1, PHI+1.0).normalize());
	GDFVectors.push(new THREE.Vector3(0, -1, PHI+1.0).normalize());
	GDFVectors.push(new THREE.Vector3(PHI+1.0, 0, 1).normalize());
	GDFVectors.push(new THREE.Vector3(-PHI-1.0, 0, 1).normalize());
	GDFVectors.push(new THREE.Vector3(1, PHI+1.0, 0).normalize());
	GDFVectors.push(new THREE.Vector3(-1, PHI+1.0, 0).normalize());
	GDFVectors.push(new THREE.Vector3(0, PHI, 1).normalize());
	GDFVectors.push(new THREE.Vector3(0, -PHI, 1).normalize());
	GDFVectors.push(new THREE.Vector3(1, 0, PHI).normalize());
	GDFVectors.push(new THREE.Vector3(-1, 0, PHI).normalize());
	GDFVectors.push(new THREE.Vector3(PHI, 1, 0).normalize());
	GDFVectors.push(new THREE.Vector3(-PHI, 1, 0).normalize());

	var GDFVectorsFlattened = [];
	for (var n=0; n<GDFVectors.length; n++) 
	{
		var v = GDFVectors[n];
		GDFVectorsFlattened.push(v.x);
		GDFVectorsFlattened.push(v.y);
		GDFVectorsFlattened.push(v.z);
	}

	traceProgram.uniform3Fv("GDFVectors", GDFVectorsFlattened);

	traceProgram.uniformF("_size", this._settings.size);
	traceProgram.uniformF("_bulge", this._settings.bulge);
	traceProgram.uniformF("_spacing", this._settings.spacing);
	traceProgram.uniformF("_width",  this._settings.width);
	traceProgram.uniformF("_height", this._settings.height);
	traceProgram.uniformF("_depth",  this._settings.depth);
	traceProgram.uniformF("_offset",  this._settings.offset);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
LatticeScene.prototype.getScale = function()
{
	return 1.0;
}


/*
LatticeScene.prototype.getBox = function()
{
	var min = new THREE.Vector3(-100, -100, -100);
	var max = new THREE.Vector3(100, 100, 100);
	return new THREE.Box3(min, max);
}
*/


// Initial laser position and direction defaults for this scene
LatticeScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(0.00000, 0.00000, 17.7086));
	laser.setDirection(new THREE.Vector3(2.22045e-16, 0.120146, -0.992756));
	laser.setEmissionRadius(0.00000);
	laser.setEmissionSpreadAngle(2.00000);
	controls.target.set(-1.17646, -0.597977, 1.94572);
	camera.position.set(21.2816, 5.70411, 27.1872);
}


// set up gui and callbacks for this scene
LatticeScene.prototype.initGui = function(parentFolder)
{
	this.typeItem = parentFolder.add(this._settings, 'shape', this._settings.shapes);
	this.typeItem.onChange( function(value) { snelly.reset(); } );

	this.rItem = parentFolder.add(this._settings, 'size', 0.1, 10.0);
	this.bItem = parentFolder.add(this._settings, 'bulge', 0.0, 1.0);
	this.nItem = parentFolder.add(this._settings, 'spacing', 0.1, 10.0);
	this.wItem = parentFolder.add(this._settings, 'width', 1.0, 100.0);
	this.hItem = parentFolder.add(this._settings, 'height', 1.0, 100.0);
	this.dItem = parentFolder.add(this._settings, 'depth', 1.0, 100.0);
	this.oItem = parentFolder.add(this._settings, 'offset', 0.0, 1.0);
	
	this.rItem.onChange( function(value) { snelly.reset(); } );
	this.bItem.onChange( function(value) { snelly.reset(); } );
	this.nItem.onChange( function(value) { snelly.reset(); } );
	this.wItem.onChange( function(value) { snelly.reset(); } );
	this.hItem.onChange( function(value) { snelly.reset(); } );
	this.dItem.onChange( function(value) { snelly.reset(); } );
	this.oItem.onChange( function(value) { snelly.reset(); } );
}


LatticeScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.typeItem);
	parentFolder.remove(this.rItem);
	parentFolder.remove(this.bItem);
	parentFolder.remove(this.nItem);
	parentFolder.remove(this.wItem);
	parentFolder.remove(this.hItem);
	parentFolder.remove(this.dItem);
	parentFolder.remove(this.oItem);
}











