

////////////////////////////////////////////////////////
// Material
////////////////////////////////////////////////////////

function Material(name, desc)
{
	this._name = name;
	this._desc = desc;
}

Material.prototype.getName = function()
{
	return this._name;
}

Material.prototype.getDesc = function()
{
	return this._desc;
}


////////////////////////////////////////////////////////
// Metals
////////////////////////////////////////////////////////

function Metal(name, desc)
{
	Material.call(this, name, desc);
}

Metal.prototype = Object.create(Material.prototype);

Metal.prototype.sample = function()
{
	return `
				float SAMPLE(inout vec3 X, inout vec3 D, vec3 N, float lnm, inout vec4 rnd)
				{                                                          
					return sampleMetal(X, D, N, IOR(lnm), K(lnm), rnd);    
				}
	`;
}


// Linear Metal
//    An extremely crude model -- we just draw a line between the values at 
//    of n and k at 400nm and 700nm.
function LinearMetal(name, desc, n400, k400, n700, k700)
{
	Metal.call(this, name, desc);
	this.n400 = n400;
	this.n700 = n700;
	this.k400 = k400;
	this.k700 = k700;
}

LinearMetal.prototype = Object.create(Metal.prototype);

LinearMetal.prototype.ior = function()
{
	// Defines GLSL functions which take wavelength (in nanometres) and return ior and k
	return `
				uniform float _n400;
				uniform float _n700;
				uniform float _k400;
				uniform float _k700;
				float IOR(float lnm)
				{                                                       
					return abs(_n400 + (lnm-400.0)*(_n700-_n400)/300.0); 
				}                                                       
				float K(float lnm)                                      
				{                                                       
					return abs(_k400 + (lnm-400.0)*(_k700-_k400)/300.0); 
				}
	`;
}

LinearMetal.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_n400", this.n400);
	traceProgram.uniformF("_n700", this.n700);
	traceProgram.uniformF("_k400", this.k400);
	traceProgram.uniformF("_k400", this.k700);
}

// set up gui and callbacks for this material
LinearMetal.prototype.initGui = function(parentFolder)
{

}

LinearMetal.prototype.eraseGui = function(parentFolder)
{

}


////////////////////////////////////////////////////////
// Dielectrics
////////////////////////////////////////////////////////

function Dielectric(name, desc)
{
	Material.call(this, name, desc);
}

Dielectric.prototype = Object.create(Material.prototype);

Dielectric.prototype.sample = function()
{
	return `
				float SAMPLE(inout vec3 X, inout vec3 D, vec3 N, float lnm, inout vec4 rnd)
				{                                                          
					return sampleDielectric(X, D, N, IOR(lnm), rnd);       
				}
	`;
}

//
// Simplest (but unphysical) model with no wavelength dependence
//
function ConstantDielectric(name, desc, iorVal) 
{
	Dielectric.call(this, name, desc);
	this.iorVal = iorVal;
}

ConstantDielectric.prototype = Object.create(Dielectric.prototype);

ConstantDielectric.prototype.ior = function()
{
	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	return `
				uniform float _iorVal;
				float IOR(float lnm)  
				{                     
					return _iorVal;   
				}
	`;
}

ConstantDielectric.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_iorVal", this.iorVal);
}

// set up gui and callbacks for this material
ConstantDielectric.prototype.initGui = function(parentFolder)
{
	this.iorItem = parentFolder.add(this, 'iorVal', 0.0, 5.0);
	this.iorItem.onChange( function(value) { renderer.reset(); } );
}

ConstantDielectric.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.iorItem);
}


//
// The standard Sellmeier model for dielectrics.
//
function SellmeierDielectric(name, desc, C1, C2, C3, C4, C5, C6, C7) 
{
	Dielectric.call(this, name, desc);
	this.C1 = C1;
	this.C2 = C2;
	this.C3 = C3;
	this.C4 = C4;
	this.C5 = C5;
	this.C6 = C6;
	this.C7 = C7;
}

SellmeierDielectric.prototype = Object.create(Dielectric.prototype);

SellmeierDielectric.prototype.ior = function()
{
	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	return `
				uniform float _C1;    
				uniform float _C2;    
				uniform float _C3;    
				uniform float _C4;    
				uniform float _C5;    
				uniform float _C6;    
				uniform float _C7;    
				float IOR(float lnm) 
				{                                                                                            
					float lmum = 1.0e-3*lnm;                                                                      
					float l2 = lmum*lmum;                                                                               
					float n2 = 1.0 + _C1 + _C2*l2/(l2 - _C3*_C3) + _C4*l2/(l2 - _C5*_C5) + _C6*l2/(l2 - _C7*_C7); 
					return sqrt(abs(n2));                                                                     
				}
	`;
}

SellmeierDielectric.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_C1", this.C1);
	traceProgram.uniformF("_C2", this.C2);
	traceProgram.uniformF("_C3", this.C3);
	traceProgram.uniformF("_C4", this.C4);
	traceProgram.uniformF("_C5", this.C5);
	traceProgram.uniformF("_C6", this.C6);
	traceProgram.uniformF("_C7", this.C7);
}

// set up gui and callbacks for this material
SellmeierDielectric.prototype.initGui = function(parentFolder)
{

}


SellmeierDielectric.prototype.eraseGui = function(parentFolder)
{

}


////////////////////////////////////////////////////////////////////////
// Instantiate a variety of material models
////////////////////////////////////////////////////////////////////////

function createMaterials(renderer)
{
	// Dielectrics
	renderer.addMaterial( new ConstantDielectric("constant", "Constant IOR dielectric", 1.5) ); 
	renderer.addMaterial( new SellmeierDielectric("diamond", "Diamond",       0.0, 0.3306,     0.175,         4.3356,      0.1060,       0.0,        0.0       ) );
	
	var glass = new SellmeierDielectric("glass_bk7", "Glass (BK7)", 0.0, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653);
	renderer.addMaterial(glass);

	// Metals
	renderer.addMaterial( new LinearMetal("aluminium", "Aluminium",  0.46555, 4.7121, 1.6620, 8.0439) );
	renderer.addMaterial( new LinearMetal("gold", "Gold",            1.5275, 1.8394, 0.16918, 3.8816) );



}


