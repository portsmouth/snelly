

////////////////////////////////////////////////////////
// Material
////////////////////////////////////////////////////////

function Material(name, desc)
{
	this._name = name;
	this._desc = desc;
	this.roughness = 0.001;
}

Material.prototype.getName = function()
{
	return this._name;
}

Material.prototype.getDesc = function()
{
	return this._desc;
}

Material.prototype.getRoughness = function()
{
	return this.roughness;
}

Material.prototype.initGui  = function(parentFolder) 
{ 
	this.roughnessItem = parentFolder.add(this, 'roughness', 0.0, 1.0);
	this.roughnessItem.onChange( function(value) { snelly.reset(); } );
}

Material.prototype.eraseGui = function(parentFolder) 
{ 
	parentFolder.remove(this.roughnessItem);
}

Material.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("roughness", this.roughness);
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
				float SAMPLE(inout vec3 X, inout vec3 D, vec3 N, float wavelength_nm, inout vec4 rnd)
				{                                                          
					return sampleMetal(X, D, N, IOR(wavelength_nm), K(wavelength_nm), rnd);    
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
				float IOR(float wavelength_nm)
				{                                                       
					return abs(_n400 + (wavelength_nm-400.0)*(_n700-_n400)/300.0); 
				}                                                       
				float K(float wavelength_nm)                                      
				{                                                       
					return abs(_k400 + (wavelength_nm-400.0)*(_k700-_k400)/300.0); 
				}
	`;
}

LinearMetal.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_n400", this.n400);
	traceProgram.uniformF("_n700", this.n700);
	traceProgram.uniformF("_k400", this.k400);
	traceProgram.uniformF("_k400", this.k700);
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
LinearMetal.prototype.initGui = function(parentFolder)
{

}

LinearMetal.prototype.eraseGui = function(parentFolder)
{

}

LinearMetal.prototype.initGui  = function(parentFolder) { Material.prototype.initGui.call(this, parentFolder) }
LinearMetal.prototype.eraseGui = function(parentFolder) { Material.prototype.eraseGui.call(this, parentFolder) }


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
				float SAMPLE(inout vec3 X, inout vec3 D, vec3 N, float wavelength_nm, inout vec4 rnd)
				{                                                          
					return sampleDielectric(X, D, N, IOR(wavelength_nm), rnd);       
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
				float IOR(float wavelength_nm)  
				{                     
					return _iorVal;   
				}
	`;
}

ConstantDielectric.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF("_iorVal", this.iorVal);
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
ConstantDielectric.prototype.initGui = function(parentFolder)
{
	this.iorItem = parentFolder.add(this, 'iorVal', 0.0, 5.0);
	this.iorItem.onChange( function(value) { snelly.reset(); } );
	Material.prototype.initGui.call(this, parentFolder);
}

ConstantDielectric.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.iorItem);
	Material.prototype.eraseGui.call(this, parentFolder)
}


// The standard Sellmeier model for dielectrics (model 1 at refractiveindex.info)
function SellmeierDielectric(name, desc, coeffs) 
{
	Dielectric.call(this, name, desc);
	this.coeffs = coeffs;
}

SellmeierDielectric.prototype = Object.create(Dielectric.prototype);

SellmeierDielectric.prototype.ior = function()
{
	var numTerms = (this.coeffs.length - 1)/2;
	var IOR_FORMULA = `1.0 + _C1 `;
	for (var t=1; t<=numTerms; ++t)
	{
		IOR_FORMULA += `+ _C${2*t}*l2/(l2 - _C${2*t+1}*_C${2*t+1})`;
	}

	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	var uniforms = '';
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		uniforms += `uniform float _C${n};\n`
	}
	var code = `${uniforms}    
	float IOR(float wavelength_nm) 
	{                                                                                            
		float wavelength_um = 1.0e-3*wavelength_nm;                                                                      
		float l2 = wavelength_um*wavelength_um;                                                                               
		float n2 = ${IOR_FORMULA}; 
		return max(sqrt(abs(n2)), 1.0e-3);                                                                     
	}`;

	return code;
}

SellmeierDielectric.prototype.syncShader = function(traceProgram)
{
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		traceProgram.uniformF(`_C${n}`, this.coeffs[n-1]);
	}
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
SellmeierDielectric.prototype.initGui  = function(parentFolder) { Material.prototype.initGui.call(this, parentFolder) }
SellmeierDielectric.prototype.eraseGui = function(parentFolder) { Material.prototype.eraseGui.call(this, parentFolder) }


// The standard Sellmeier model for dielectrics (model 2 at refractiveindex.info)
// coeffs array must have an odd number of elements (the constant, plus a pair per 'pole' term)
function Sellmeier2Dielectric(name, desc, coeffs) 
{
	Dielectric.call(this, name, desc);
	this.coeffs = coeffs;
}

Sellmeier2Dielectric.prototype = Object.create(Dielectric.prototype);

Sellmeier2Dielectric.prototype.ior = function()
{
	var numTerms = (this.coeffs.length - 1)/2;
	var IOR_FORMULA = `1.0 + _C1 `;
	for (var t=1; t<=numTerms; ++t)
	{
		IOR_FORMULA += `+ _C${2*t}*l2/(l2 - _C${2*t+1})`;
	}

	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	var uniforms = '';
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		uniforms += `uniform float _C${n};\n`
	}
	var code = `${uniforms}    
	float IOR(float wavelength_nm) 
	{                                                                                            
		float wavelength_um = 1.0e-3*wavelength_nm;                                                                      
		float l2 = wavelength_um*wavelength_um;                                                                               
		float n2 = ${IOR_FORMULA}; 
		return max(sqrt(abs(n2)), 1.0e-3);                                                                     
	}`;

	return code;
}

Sellmeier2Dielectric.prototype.syncShader = function(traceProgram)
{
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		traceProgram.uniformF(`_C${n}`, this.coeffs[n-1]);
	}
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
Sellmeier2Dielectric.prototype.initGui  = function(parentFolder) { Material.prototype.initGui.call(this, parentFolder) }
Sellmeier2Dielectric.prototype.eraseGui = function(parentFolder) { Material.prototype.eraseGui.call(this, parentFolder) }



// Model 4 at Polyanskiy's refractiveindex.info:
function PolyanskiyDielectric(name, desc, coeffs) 
{
	Dielectric.call(this, name, desc);
	this.C1 = coeffs[0];
	this.C2 = coeffs[1];
	this.C3 = coeffs[2];
	this.C4 = coeffs[3];
	this.C5 = coeffs[4];
}

PolyanskiyDielectric.prototype = Object.create(Dielectric.prototype);

PolyanskiyDielectric.prototype.ior = function()
{
	var IOR_FORMULA = ' _C1 + _C2*pow(l, _C3)/(l*l - pow(_C4, _C5))'; 

	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	var code = `
	uniform float _C1;
	uniform float _C2;
	uniform float _C3;
	uniform float _C4;
	uniform float _C5;
	float IOR(float wavelength_nm) 
	{                                                                                            
		float wavelength_um = 1.0e-3*wavelength_nm;                                                                      
		float l = wavelength_um;                                                                               
		float n2 = ${IOR_FORMULA}; 
		return max(sqrt(abs(n2)), 1.0e-3);                                                                     
	}`;

	return code;
}

PolyanskiyDielectric.prototype.syncShader = function(traceProgram)
{
	traceProgram.uniformF('_C1', this.C1);
	traceProgram.uniformF('_C2', this.C2);
	traceProgram.uniformF('_C3', this.C3);
	traceProgram.uniformF('_C4', this.C4);
	traceProgram.uniformF('_C5', this.C5);
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
PolyanskiyDielectric.prototype.initGui  = function(parentFolder) { Material.prototype.initGui.call(this, parentFolder) }
PolyanskiyDielectric.prototype.eraseGui = function(parentFolder) { Material.prototype.eraseGui.call(this, parentFolder) }


// Cauchy model for dielectrics (model 5 at refractiveindex.info)
function CauchyDielectric(name, desc, coeffs) 
{
	Dielectric.call(this, name, desc);
	this.coeffs = coeffs;
}

CauchyDielectric.prototype = Object.create(Dielectric.prototype);

CauchyDielectric.prototype.ior = function()
{
	var numTerms = (this.coeffs.length - 1)/2;
	var IOR_FORMULA = `_C1`;
	for (var t=1; t<=numTerms; ++t)
	{
		IOR_FORMULA += ` + _C${2*t}*pow(l, _C${2*t+1})`;
	}

	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	var uniforms = '';
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		uniforms += `uniform float _C${n};\n`
	}
	var code = `${uniforms}    
	float IOR(float wavelength_nm) 
	{                                                                                            
		float wavelength_um = 1.0e-3*wavelength_nm;                                                                      
		float l = wavelength_um;                                                                               
		float n = ${IOR_FORMULA}; 
		return max(n, 1.0e-3);                                                                     
	}`;

	return code;
}

CauchyDielectric.prototype.syncShader = function(traceProgram)
{
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		traceProgram.uniformF(`_C${n}`, this.coeffs[n-1]);
	}
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
CauchyDielectric.prototype.initGui  = function(parentFolder) { Material.prototype.initGui.call(this, parentFolder) }
CauchyDielectric.prototype.eraseGui = function(parentFolder) { Material.prototype.eraseGui.call(this, parentFolder) }



// Gases (model 6 at refractiveindex.info)
function Gas(name, desc, coeffs) 
{
	Dielectric.call(this, name, desc);
	this.coeffs = coeffs;
}

Gas.prototype = Object.create(Dielectric.prototype);

Gas.prototype.ior = function()
{
	var numTerms = (this.coeffs.length - 1)/2;
	var IOR_FORMULA = `1.0 + _C1 `;
	for (var t=1; t<=numTerms; ++t)
	{
		IOR_FORMULA += `+ _C${2*t}/(_C${2*t+1} - invl2)`;
	}

	// Defines a GLSL function which takes wavelength (in micrometres) and returns ior
	var uniforms = '';
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		uniforms += `uniform float _C${n};\n`
	}
	var code = `${uniforms}    
	float IOR(float wavelength_nm) 
	{                                                                                            
		float wavelength_um = 1.0e-3*wavelength_nm;                                                                      
		float invl2 = 1.0 / (wavelength_um*wavelength_um);                                                                               
		float n2 = ${IOR_FORMULA}; 
		return max(sqrt(abs(n2)), 1.0e-3);                                                                     
	}`;

	return code;
}

Gas.prototype.syncShader = function(traceProgram)
{
	for (var n=1; n<=this.coeffs.length; ++n)
	{
		traceProgram.uniformF(`_C${n}`, this.coeffs[n-1]);
	}
	Material.prototype.syncShader.call(this, traceProgram);
}

// set up gui and callbacks for this material
Gas.prototype.initGui  = function(parentFolder) { Material.prototype.initGui.call(this, parentFolder) }
Gas.prototype.eraseGui = function(parentFolder) { Material.prototype.eraseGui.call(this, parentFolder) }


