

var Snelly = function(sceneObj)
{
	this.initialized = false; 
	snelly = this;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;
	window.addEventListener( 'resize', this, false );
	this.container = document.getElementById('container');

	// Text canvas
	var text_canvas = document.getElementById("text-canvas");
	this.textCtx = text_canvas.getContext("2d");
	this.onSnellyLink = false;

	// Setup THREE.js GL camera and orbit controls for it
	var VIEW_ANGLE = 45;
	var ASPECT = this.width / this.height ;
	var NEAR = 0.05;
	var FAR = 1000;
	this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
	this.controls = new THREE.OrbitControls(this.camera, this.container);
	this.controls.zoomSpeed = 2.0;
	this.controls.addEventListener( 'change', camChanged );

	this.gui = null;

	// Load scene
	this.loadScene(sceneObj);

	// Instantiate materials
	this.dielectrics = {}
	this.metals = {}
	this.dielectricObj = null;
	this.metalObj = null;
	{
		// Dielectrics
		this.addDielectric( new ConstantDielectric("Constant IOR dielectric", "", 1.5) ); 
		this.addDielectric( new SellmeierDielectric("Glass (BK7)", "",       [0.0, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945,  103.560653]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (K7)", "",       [0.0, 1.1273555,  0.00720341707, 0.124412303, 0.0269835916, 0.827100531, 100.384588]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (F5)", "",       [0.0, 1.3104463,  0.00958633048, 0.19603426,  0.0457627627, 0.96612977,  115.011883]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (LAFN7)", "",    [0.0, 1.66842615, 0.0103159999,  0.298512803, 0.0469216348, 1.0774376,   82.5078509]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (LASF35)", "",   [0.0, 2.45505861, 0.0135670404,  0.453006077, 0.054580302,  2.3851308,   167.904715]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (N-LAK33A)", "", [0.0, 1.44116999, 0.00680933877, 0.571749501, 0.0222291824, 1.16605226,  80.9379555]) );
		this.addDielectric( new SellmeierDielectric("Glass (N-FK51A)", "",   [0.0, 0.97124781, 0.00472301995, 0.216901417, 0.0153575612, 0.90465166,  168.68133]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (SF4)", "",      [0.0, 1.61957826, 0.0125502104,  0.339493189, 0.0544559822, 1.02566931,  117.652222]) );
		this.addDielectric( new Sellmeier2Dielectric("Glass (SF67)", "",     [0.0, 1.97464225, 0.0145772324,  0.467095921, 0.0669790359, 2.43154209,  157.444895]) );
		this.addDielectric( new Sellmeier2Dielectric("Water", "",            [0.0,        5.67252e-1, 5.08555046e-3, 1.736581e-1, 1.8149386e-2, 2.12153e-2, 2.61726e-2, 1.1384932e-1, 1.073888e1]) );
		this.addDielectric( new Sellmeier2Dielectric("Ethanol", "",          [0.0,        0.83189,    0.00930,       -0.15582,    -49.45200]) );
		this.addDielectric( new Sellmeier2Dielectric("Polycarbonate", "",    [0.0,        0.83189,    0.00930,       -0.15582,    -49.45200]) );
		this.addDielectric( new CauchyDielectric("Glycerol", "",             [1.45797, 0.00598, -2, -0.00036, -4]) );
		this.addDielectric( new CauchyDielectric("Liquid Crystal (E7)", "",  [1.4990,  0.0072,  -2,  0.0003,  -4]) );
		this.addDielectric( new SellmeierDielectric("Diamond", "",           [0.0,        0.3306,     0.175,         4.3356,      0.1060]) );
		this.addDielectric( new SellmeierDielectric("Quartz", "",            [0.0, 0.6961663, 0.0684043, 0.4079426, 0.1162414, 0.8974794, 9.896161]) );
		this.addDielectric( new SellmeierDielectric("Fused Silica", "",      [0.0,        0.6961663,  0.0684043,     0.4079426,  0.1162414, 0.8974794, 9.896161]) );
		this.addDielectric( new SellmeierDielectric("Sapphire", "",          [0.0,        1.5039759,  0.0740288,     0.55069141, 0.1216529, 6.5927379, 20.072248]) );
		this.addDielectric( new SellmeierDielectric("Sodium Chloride", "",   [0.00055,    0.19800,    0.050,         0.48398,     0.100,        0.38696,   0.128]) );
		this.addDielectric( new PolyanskiyDielectric("Proustite", "",        [7.483, 0.474, 0.0, 0.09, 1.0]) );
		this.addDielectric( new PolyanskiyDielectric("Rutile (Titanium Dioxide)", "", [5.913, 0.2441, 0.0, 0.0803, 1.0]) );
		this.addDielectric( new PolyanskiyDielectric("Silver Chloride", "", [4.00804, 0.079086, 0.0, 0.04584, 1.0]) );

		// @todo: Metals. Need a much better approximation than 'LinearMetal' for plausible metals 
		this.addMetal( new LinearMetal("Aluminium", "",  0.46555, 4.7121, 1.6620, 8.0439) );
		this.addMetal( new LinearMetal("Gold", "",       1.5275, 1.8394, 0.16918, 3.8816) );
	}

	// Load the initial material
	this.loadDielectric("Glass (LASF35)");
	this.loadMetal("Gold");

	// Instantiate distance field pathtracer
	this.pathtracer = new Pathtracer();

	// Do initial resize:
	this.resize();
	
	// Instantiate spectra
	{
		// Spectrum initialization
		this.spectra = {}
		this.SPECTRUM_SAMPLES = 1024;
		this.spectrumObj = null;
		this.LAMBDA_MIN = 390.0;
	    this.LAMBDA_MAX = 750.0;
		var wToRgb = wavelengthToRgbTable();
		this.wavelengthToRgb = new GLU.Texture(wToRgb.length/4, 1, 4, true,  true, true, wToRgb);
		this.emissionIcdf    = new GLU.Texture(4*this.SPECTRUM_SAMPLES, 1, 1, true, true, true, null);

		this.addSpectrum( new FlatSpectrum("flat", "Flat spectrum", 400.0, 700.0) );
		this.addSpectrum( new BlackbodySpectrum("blackbody", "Blackbody spectrum", 6000.0) );
		this.addSpectrum( new MonochromaticSpectrum("monochromatic", "Monochromatic spectrum", 650.0) ); 

		this.loadSpectrum("blackbody");
	}

	// Create dat gui
	this.gui = new GUI();

	// Setup keypress and mouse events
	window.addEventListener('keydown', function(event) 
	{
		var charCode = (event.which) ? event.which : event.keyCode;
		switch (charCode)
		{
			case 122: // F11 key: go fullscreen
				var element	= document.body;
				if      ( 'webkitCancelFullScreen' in document ) element.webkitRequestFullScreen();
				else if ( 'mozCancelFullScreen'    in document ) element.mozRequestFullScreen();
				else console.assert(false);
				break;

			case 70: // F key: reset cam  (@todo)
				break;

			case 82: // R key: reset scene 
				// @todo
				break;

			case 72: // H key: hide dat gui
				snelly.getGUI().visible = !(snelly.getGUI().visible);
				break;

			case 67: // C key: dev tool to dump cam and laser details, for setting scene defaults
				var t = snelly.controls.target;
				var c = snelly.camera.position;
				console.log(`controls.target.set(${t.x.toPrecision(6)}, ${t.y.toPrecision(6)}, ${t.z.toPrecision(6)});`);
				console.log(`camera.position.set(${c.x.toPrecision(6)}, ${c.y.toPrecision(6)}, ${c.z.toPrecision(6)});`);
				break;
			
		}
	}, false);

	window.addEventListener( 'mousemove', this, false );
	window.addEventListener( 'mousedown', this, false );
	window.addEventListener( 'mouseup',   this, false );
	window.addEventListener( 'contextmenu',   this, false );
	window.addEventListener( 'click', this, false );

	this.initialized = true; 
}

Snelly.prototype.handleEvent = function(event)
{
	switch (event.type)
	{
		case 'resize':      this.resize();  break;
		case 'mousemove':   this.onDocumentMouseMove(event);  break;
		case 'mousedown':   this.onDocumentMouseDown(event);  break;
		case 'mouseup':     this.onDocumentMouseUp(event);    break;
		case 'contextmenu': this.onDocumentRightClick(event); break;
		case 'click':       this.onClick(event);  break;
	}
}

Snelly.prototype.getPathtracer = function()
{
	return this.pathtracer;
}

Snelly.prototype.getGUI = function()
{
	return this.gui;
}

Snelly.prototype.getCamera = function()
{
	return this.camera;
}


//
// Scene management
//

Snelly.prototype.getScene = function()
{
	return this.sceneObj;
}

Snelly.prototype.loadScene = function(sceneObj)
{
	this.sceneObj = sceneObj;
	this.sceneObj.init(this.controls, this.camera, this.laser);

	var gui = this.getGUI();
	if (gui)
	{
		gui.emissionRadiusControl.max(4.0*this.sceneObj.getScale());
	}
	
	// Camera frustum update
	this.camera.near = Math.max(1.0e-4, 1.0e-2*this.sceneObj.getScale());
	this.camera.far  = Math.max(1.0e4,   1.0e4*this.sceneObj.getScale());
	this.controls.update();	
	this.reset();
}


//
// Spectrum management
//

// emission spectrum management
Snelly.prototype.addSpectrum = function(spectrumObj)
{
	this.spectra[spectrumObj.getName()] = spectrumObj;
}

Snelly.prototype.getSpectra = function()
{
	return this.spectra;
}

Snelly.prototype.loadSpectrum = function(spectrumName)
{
	this.spectrumObj = this.spectra[spectrumName];
	var inverseCDF = this.spectrumObj.inverseCDF(this.LAMBDA_MIN, this.LAMBDA_MAX, this.SPECTRUM_SAMPLES);
	this.emissionIcdf.bind(0);
    this.emissionIcdf.copy(inverseCDF);
    this.reset();
}

Snelly.prototype.getLoadedSpectrum = function()
{
	return this.spectrumObj;
}

//
// Material management
//
Snelly.prototype.addDielectric = function(materialObj)
{
	this.dielectrics[materialObj.getName()] = materialObj;
}

Snelly.prototype.getDielectrics = function()
{
	return this.dielectrics;
}

Snelly.prototype.addMetal = function(materialObj)
{
	this.metals[materialObj.getName()] = materialObj;
}

Snelly.prototype.getMetals = function()
{
	return this.metals;
}


Snelly.prototype.loadDielectric = function(dielectricName)
{
	this.dielectricObj = this.dielectrics[dielectricName];
	this.reset();
}

Snelly.prototype.loadMetal = function(metalName)
{
	this.metalObj = this.metals[metalName];
	this.reset();
}

Snelly.prototype.getLoadedDielectric = function()
{
	return this.dielectricObj;
}

Snelly.prototype.getLoadedMetal = function()
{
	return this.metalObj;
}


// Renderer reset on camera update
Snelly.prototype.reset = function(no_recompile = false)
{	
	if (!this.initialized) return;

	this.pathtracer.reset(no_recompile);

	this.gui.sync();

	this.render();
}

// Render all 
Snelly.prototype.render = function()
{
	if (!this.initialized) return;
	if (this.sceneObj == null) return;
	var gl = GLU.gl;
	gl.viewport(0, 0, this.width, this.height);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.disable(gl.DEPTH_TEST);

	// Render distance field surface via pathtracing
	this.pathtracer.render();

	// Update HUD text canvas	
	this.textCtx.textAlign = "left";   	// This determines the alignment of text, e.g. left, center, right
	this.textCtx.textBaseline = "middle";	// This determines the baseline of the text, e.g. top, middle, bottom
	this.textCtx.font = '12px monospace';	// This determines the size of the text and the font family used
	this.textCtx.clearRect(0, 0, this.textCtx.canvas.width, this.textCtx.canvas.height);
	this.textCtx.globalAlpha = 0.95;
	if (snelly.getGUI().visible)
	{
	  	if (this.onSnellyLink) this.textCtx.fillStyle = "#ff5500";
	  	else                   this.textCtx.fillStyle = "#ffff00";
	  	this.textCtx.fillText('Snelly renderer', 14, 20);
	  	this.textCtx.fillStyle = "#aaaaff";
	  	//this.textCtx.fillText('ray count:    ' + (lsStats.rayCount/1.0e6).toPrecision(3) + 'M', 14, 35);
	  	//this.textCtx.fillText('waves traced: ' + lsStats.wavesTraced,    14, 50);
	}
}

Snelly.prototype.resize = function()
{
	var width = window.innerWidth;
	var height = window.innerHeight;
	this.width = width;
	this.height = height;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = width;
	render_canvas.height = height;

	var text_canvas = document.getElementById("text-canvas");
	text_canvas.width  = width;
	text_canvas.height = height

	this.camera.aspect = width / height;
	this.camera.updateProjectionMatrix();

	this.pathtracer.resize(width, height);

	if (this.initialized)
		this.render();
}


Snelly.prototype.onClick = function(event)
{
	if (this.onSnellyLink) 
	{
    	window.location = "https://github.com/portsmouth/snelly";
    }
	event.preventDefault();
}

Snelly.prototype.onDocumentMouseMove = function(event)
{
	// Check whether user is trying to click the Snelly home link
	var textCtx = this.textCtx;
	var x = event.pageX;
    var y = event.pageY;
	var linkWidth = this.textCtx.measureText('Snelly renderer').width;
	if (x>14 && x<14+linkWidth &&
		y>15 && y<25 )
	{
		this.onSnellyLink = true;
	}
	else
	{
		this.onSnellyLink = false;
	}

	console.log('mouse move');
	this.controls.update();
	event.preventDefault();
}

Snelly.prototype.onDocumentMouseDown = function(event)
{
	this.controls.update();
	//event.preventDefault();
}

Snelly.prototype.onDocumentMouseUp = function(event)
{
	this.controls.update();
	//event.preventDefault();
}

Snelly.prototype.onDocumentRightClick = function(event)
{
	this.controls.update();
	event.preventDefault();
	if (event.altKey) return; // don't pick if alt-right-clicking (panning)

	var xPick =   (( event.clientX - window.offsetLeft ) / window.width)*2 - 1;
	var yPick = - (( event.clientY - window.offsetTop ) / window.height)*2 + 1;

	var pickedPoint = this.pathtracer.pick(xPick, yPick);
	if (pickedPoint == null)
	{
		// unset?
		return;
	} 

	var no_recompile = true;
	this.reset(no_recompile);
}

function camChanged()
{
	var no_recompile = true;
	snelly.reset(no_recompile);
	snelly.render();
}





