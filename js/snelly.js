

var Snelly = function(sceneObj)
{
	this.initialized = false; 
	snelly = this;

	// @todo: should create the required canvas elements here
	//        programmatically, rather than relying on them existing in the document.

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

	// Setup THREE.js orbit camera
	var VIEW_ANGLE = 45; // @todo: fov should be under user control
	var ASPECT = this.width / this.height;
	var NEAR = 0.05;
	var FAR = 1000;
	this.glCamera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
	this.camControls = new THREE.OrbitControls(this.glCamera, this.container);
	this.camControls.zoomSpeed = 2.0;
	this.camControls.addEventListener( 'change', camChanged );
	this.camControls.target.set(0.0, 0.0, 0.0);
	this.glCamera.position.set(1.0, 1.0, 1.0);

	this.gui = null;

	// Load scene
	this.loadScene(sceneObj);

	// Instantiate materials
	this.materials = new Materials();

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
		var wToXYZ = wavelengthToXYZTable();
		this.wavelengthToXYZ = new GLU.Texture(wToXYZ.length/4, 1, 4, true,  true, true, wToXYZ);
		this.emissionIcdf    = new GLU.Texture(4*this.SPECTRUM_SAMPLES, 1, 1, true, true, true, null);

		this.addSpectrum( new FlatSpectrum("flat", "Flat spectrum", 400.0, 700.0) );
		this.addSpectrum( new BlackbodySpectrum("blackbody", "Blackbody spectrum", 6000.0) );
		this.addSpectrum( new MonochromaticSpectrum("monochromatic", "Monochromatic spectrum", 650.0) ); 

		this.loadSpectrum("blackbody");
	}

	// Create dat gui
	this.gui = new GUI();

	// Setup keypress and mouse events
	window.addEventListener( 'mousemove', this, false );
	window.addEventListener( 'mousedown', this, false );
	window.addEventListener( 'mouseup',   this, false );
	window.addEventListener( 'contextmenu',   this, false );
	window.addEventListener( 'click', this, false );
	window.addEventListener( 'keydown', this, false );

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
		case 'keydown':     this.onKeydown(event);  break;
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
	return this.glCamera;
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
	this.sceneObj.init(this.camControls, this.glCamera);

	// @todo: if no init function provided, init cam pos 
	//        according to scene scale
	//camControls.target.set(0.0, 0.0, 0.0);
	//camera.position.set(1.0, 1.0, 1.0);
	
	// Camera frustum update
	this.glCamera.near = Math.max(1.0e-4, 1.0e-2*this.sceneObj.getScale());
	this.glCamera.far  = Math.max(1.0e4,   1.0e4*this.sceneObj.getScale());
	this.camControls.update();	
	this.reset();
}


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
    this.reset(true);
}

Snelly.prototype.getLoadedSpectrum = function()
{
	return this.spectrumObj;
}


// Material management
Snelly.prototype.getDielectrics = function()
{
	return this.materials.getDielectrics();
}

Snelly.prototype.getMetals = function()
{
	return this.materials.getMetals();
}

Snelly.prototype.loadDielectric = function(dielectricName)
{
	this.materials.loadDielectric(dielectricName);
	this.reset();
}

Snelly.prototype.loadMetal = function(metalName)
{
	this.materials.loadMetal(metalName);
	this.reset();
}

Snelly.prototype.getLoadedDielectric = function()
{
	return this.materials.getLoadedDielectric();
}

Snelly.prototype.getLoadedMetal = function()
{
	return this.materials.getLoadedMetal();
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
	var text_canvas = document.getElementById("text-canvas");
	render_canvas.width  = width;
	render_canvas.height = height;
	text_canvas.width  = width;
	text_canvas.height = height

	this.glCamera.aspect = width / height;
	this.glCamera.updateProjectionMatrix();

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

	this.camControls.update();
	//event.preventDefault();
}

Snelly.prototype.onDocumentMouseDown = function(event)
{
	this.camControls.update();
	//event.preventDefault();
}

Snelly.prototype.onDocumentMouseUp = function(event)
{
	this.camControls.update();
	//event.preventDefault();
}

Snelly.prototype.onDocumentRightClick = function(event)
{
	/*
	this.controls.update();
	event.preventDefault();
	if (event.altKey) return; // don't pick if alt-right-clicking (panning)
	var xPick =  (( event.clientX - window.offsetLeft ) / window.width)*2 - 1;
	var yPick = -(( event.clientY - window.offsetTop ) / window.height)*2 + 1;
	var pickedPoint = this.pathtracer.pick(xPick, yPick);
	if (pickedPoint == null)
	{
		// unset?
		return;
	} 
	var no_recompile = true;
	this.reset(no_recompile);
	*/
}

Snelly.prototype.onKeydown = function(event)
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
			this.sceneObj.init(this.camControls, this.glCamera);
			this.reset(true);
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
		
		case 87: console.log('w pressed'); break;
		case 65: console.log('a pressed'); break;
		case 83: console.log('s pressed'); break;
		case 68: console.log('d pressed'); break;
	}
}

function camChanged()
{
	var no_recompile = true;
	snelly.reset(no_recompile);
	snelly.render();
}





