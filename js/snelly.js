
/** 
* Snelly is the global object providing access to all functionality in the system.
* @constructor 
* @param {Scene} sceneObj - The user-defined scene
*/
var Snelly = function(sceneObj)
{
	this.initialized = false;
	this.rendering = false;
	this.sceneObj = sceneObj;
	snelly = this;

	let container = document.getElementById("container");
	this.container = container;

	var render_canvas = document.getElementById('render-canvas');
	this.render_canvas = render_canvas;
	this.width = render_canvas.width;
	this.height = render_canvas.height;

	var text_canvas = document.getElementById('text-canvas');
	this.text_canvas = text_canvas;
	this.textCtx = text_canvas.getContext("2d");
	this.onSnellyLink = false;
	this.onUserLink = false;

	window.addEventListener( 'resize', this, false );

	// Setup THREE.js orbit camera
	var VIEW_ANGLE = 45; // @todo: fov should be under user control
	var ASPECT = this.width / this.height;
	var NEAR = 0.05;
	var FAR = 1000;
	this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
	this.camera.lookAt(new THREE.Vector3(0.0, 0.0, 0.0));
	this.camera.position.set(1.0, 1.0, 1.0);
	this.camControls = new THREE.OrbitControls(this.camera, this.container);
	this.camControls.zoomSpeed = 2.0;
	this.camControls.addEventListener('change', camChanged);
	this.camControls.keyPanSpeed = 100.0;

	this.gui = null;
	this.guiVisible = true;

	// Instantiate materials
	this.materials = new Materials();

	// Instantiate distance field pathtracer
	this.pathtracer = new Renderer();
	this.auto_resize = true;
		
	// Spectrum initialization
	this.spectra = {}
	this.SPECTRUM_SAMPLES = 1024;
	this.spectrumObj = null;
	this.LAMBDA_MIN = 390.0;
    this.LAMBDA_MAX = 750.0;
	var wToXYZ = wavelengthToXYZTable();
	this.wavelengthToXYZ = new GLU.Texture(wToXYZ.length/4, 1, 4, true,  true, true, wToXYZ);
	this.emissionIcdf    = new GLU.Texture(4*this.SPECTRUM_SAMPLES, 1, 1, true, true, true, null);

	// Allow user to programmatically initialize the camera, materials, and renderer
	this.initScene();

	// Do initial resize:
	this.resize();

	// Create dat gui
	this.gui = new GUI(this.guiVisible);

	// Setup keypress and mouse events
	window.addEventListener( 'mousemove', this, false );
	window.addEventListener( 'mousedown', this, false );
	window.addEventListener( 'mouseup',   this, false );
	window.addEventListener( 'contextmenu',   this, false );
	window.addEventListener( 'click', this, false );
	window.addEventListener( 'keydown', this, false );

	this.initialized = true; 
}

/**
* Returns the current version number of the snelly system, in the format [1, 2, 3]
*  @returns {Array}
*/
Snelly.prototype.getVersion = function()
{
	return [1, 0, 0];
}

Snelly.prototype.handleEvent = function(event)
{
	switch (event.type)
	{
		case 'resize':      this.resize();  break;
		case 'mousemove':   this.onDocumentMouseMove(event);  break;
		case 'mousedown':   this.onDocumentMouseDown(event);  break;
		case 'mouseup':     this.onDocumentMouseUp(event);    break;
		//case 'contextmenu': this.onDocumentRightClick(event); break;
		case 'click':       this.onClick(event);  break;
		case 'keydown':     this.onKeydown(event);  break;
	}
}

/**
* Access to the Renderer object
*  @returns {Renderer} 
*/
Snelly.prototype.getRenderer = function()
{
	return this.pathtracer;
}

/**
* Access to the GUI object
*  @returns {GUI} 
*/
Snelly.prototype.getGUI = function()
{
	return this.gui;
}

/**
* Access to the camera object
* @returns {THREE.PerspectiveCamera}.
*/
Snelly.prototype.getCamera = function()
{
	return this.camera;
}

/**
* Access to the camera controller object
* @returns {THREE.OrbitControls}.
*/
Snelly.prototype.getControls = function()
{
	return this.camControls;
}

/**
* Programmatically show or hide the dat.GUI UI
* @param {Boolean} showGUI - toggle
*/
Snelly.prototype.showGUI = function(showGUI)
{
	this.guiVisible = showGUI;
}


//
// Scene management
//

Snelly.prototype.getScene = function()
{
	return this.sceneObj;
}

Snelly.prototype.initScene = function()
{
	if (typeof this.sceneObj.shader == "undefined")
	{
		GLU.fail('Scene must define a "shader" function!');
	}

	// Obtain min.max scale, if specified
	this.minScale = 1.0e-4;
	if (typeof this.sceneObj.getMinScale !== "undefined") 
	{
		this.minScale = this.sceneObj.getMinScale(); // sanity
	}
	this.minScale = Math.max(1.0e-6, this.minScale);

	this.maxScale = 1.0e2;
	if (typeof this.sceneObj.getMaxScale !== "undefined")
	{
		this.maxScale = this.sceneObj.getMaxScale();
	}
	this.maxScale = Math.min(1.0e6, this.maxScale); // sanity

	// Set initial default camera position and target based on max scale
	let po = 0.1*this.maxScale;
	this.camera.position.set(po, po, po);
	this.camControls.target.set(0.0, 0.0, 0.0);

	// Call user-defined init function	
	if (typeof this.sceneObj.init !== "undefined") 
	{
		this.sceneObj.init(this);
	}

	// Read optional scene name and URL, if provided
	this.sceneName = '';
	if (typeof this.sceneObj.getName !== "undefined") 
	{
		this.sceneName = this.sceneObj.getName();
	}

	this.sceneURL = '';
	if (typeof this.sceneObj.getURL !== "undefined") 
	{
		this.sceneURL = this.sceneObj.getURL();
	}

	// cache initial camera position to allow reset on 'F'
	this.initial_camera_position = new THREE.Vector3();
	this.initial_camera_position.copy(this.camera.position);
	this.initial_camera_target = new THREE.Vector3();
	this.initial_camera_target.copy(this.camControls.target);

	// Compile GLSL shaders
  	this.pathtracer.compileShaders();

  	// Fix renderer to width & height, if they were specified
  	if ((typeof this.pathtracer.width!=="undefined") && (typeof this.pathtracer.height!=="undefined"))
  	{
  		this.auto_resize = false;
  		this._resize(this.pathtracer.width, this.pathtracer.height);
  	}

  	// Set up blackbody spectrum sampling
  	this.spectra = [];
	this.addSpectrum( new BlackbodySpectrum("blackbody", "Blackbody spectrum", this.pathtracer.skyTemperature) );
  	this.loadSpectrum("blackbody");
	
	// Camera setup
	this.camera.near = this.minScale;
	this.camera.far  = this.maxScale;
	this.camControls.update();	
	this.reset(false);
}


Snelly.prototype.dumpScene = function()
{
	let camera = this.camera;
	let controls = this.camControls;
	let renderer = this.pathtracer;
	let materials = this.materials;

	var code = `/******* copy-pasted console output on 'O', begin *******/\n`;
	code += `
let renderer  = snelly.getRenderer();
let camera    = snelly.getCamera();
let controls  = snelly.getControls();
let materials = snelly.getMaterials();
	`;

	if (typeof this.sceneObj.initGenerator !== "undefined") 
	{
		code += this.sceneObj.initGenerator();
	}
	
	code += this.guiVisible ? `\nsnelly.showGUI(true);\n` : `\nsnelly.showGUI(false);\n`;

	code += `
/** Camera settings:
*		camera is a THREE.PerspectiveCamera object
* 		controls is a THREE.OrbitControls object
*/
camera.fov = ${camera.fov};
camera.up.set(${camera.up.x}, ${camera.up.y}, ${camera.up.z});
camera.position.set(${camera.position.x}, ${camera.position.y}, ${camera.position.z});
controls.target.set(${controls.target.x}, ${controls.target.y}, ${controls.target.z});
controls.zoomSpeed = ${controls.zoomSpeed};
controls.keyPanSpeed = ${controls.keyPanSpeed};

/** Renderer settings **/
renderer.renderMode = '${renderer.renderMode}';  // The other modes are: 'ao', 'normals'
renderer.maxBounces = ${renderer.maxBounces};
renderer.maxMarchSteps = ${renderer.maxMarchSteps};
renderer.radianceClamp = ${renderer.radianceClamp}; // (log scale)
renderer.skyPower = ${renderer.skyPower};
renderer.skyTemperature = ${renderer.skyTemperature};
renderer.exposure = ${renderer.exposure};
renderer.gamma = ${renderer.gamma};
renderer.whitepoint = ${renderer.whitepoint};
renderer.goalFPS = ${renderer.goalFPS};
renderer.minsSPPToRedraw = ${renderer.minsSPPToRedraw};

/** Material settings **/
let surface = materials.loadSurface();
surface.roughness = ${materials.loadSurface().roughness};
surface.ior = ${materials.loadSurface().ior};
surface.diffuseAlbedo = [${materials.loadSurface().diffuseAlbedo[0]}, ${materials.loadSurface().diffuseAlbedo[1]}, ${materials.loadSurface().diffuseAlbedo[2]}];
surface.specAlbedo = [${materials.loadSurface().specAlbedo[0]}, ${materials.loadSurface().specAlbedo[1]}, ${materials.loadSurface().specAlbedo[2]}];

let dielectric = materials.loadDielectric('${materials.getLoadedDielectric().getName()}');
dielectric.absorptionColor = [${materials.getLoadedDielectric().absorptionColor[0]}, ${materials.getLoadedDielectric().absorptionColor[1]}, ${materials.getLoadedDielectric().absorptionColor[2]}];
dielectric.absorptionScale = ${materials.getLoadedDielectric().absorptionScale}; // mfp in multiples of scene scale
dielectric.roughness = ${materials.getLoadedDielectric().roughness};

let metal = materials.loadMetal('${materials.getLoadedMetal().getName()}');
metal.roughness = ${materials.getLoadedMetal().roughness};

/******* copy-pasted console output on 'O', end *******/
	`;

	return code;
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

/**
* Get materials object
 * @returns {Materials}
*/
Snelly.prototype.getMaterials = function()
{
	return this.materials;
}

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

/**
* Get Surface object
 * @returns {Surface}
*/
Snelly.prototype.getSurface = function()
{
	return this.materials.loadSurface();
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
	this.rendering = true;

	if (!this.initialized) return;
	if (this.sceneObj == null) return;

	// Render distance field surface via pathtracing
	this.pathtracer.render();

	// Update HUD text canvas	
	this.textCtx.textAlign = "left";   	// This determines the alignment of text, e.g. left, center, right
	this.textCtx.textBaseline = "middle";	// This determines the baseline of the text, e.g. top, middle, bottom
	this.textCtx.font = '12px monospace';	// This determines the size of the text and the font family used
	this.textCtx.clearRect(0, 0, this.textCtx.canvas.width, this.textCtx.canvas.height);
	this.textCtx.globalAlpha = 0.95;
	if (this.guiVisible)
	{
	  	if (this.onSnellyLink) this.textCtx.fillStyle = "#ff5500";
	  	else                   this.textCtx.fillStyle = "#ffff00";
	  	let ver = this.getVersion();
	  	this.textCtx.fillText('Snelly renderer v'+ver[0]+'.'+ver[1]+'.'+ver[2], 14, 20);
	  	this.textCtx.fillStyle = "#ffccaa";
	  	this.textCtx.fillText('spp: ' + (this.pathtracer.spp).toPrecision(3), 14, 35);
	  	if (this.sceneName != '')
	  	{
	  		this.textCtx.fillStyle = "#ffaa22";
	  		this.textCtx.fillText(this.sceneName, 14, this.height-25);
	  	}
		if (this.sceneURL != '')
		{
			if (this.onUserLink) this.textCtx.fillStyle = "#aaccff";
	  		else                 this.textCtx.fillStyle = "#55aaff";
	  		this.textCtx.fillText(this.sceneURL, 14, this.height-40);
		}
	}

	this.rendering = false;
}

Snelly.prototype._resize = function(width, height)
{
	this.width = width;
	this.height = height;

	var render_canvas = this.render_canvas;
	render_canvas.width  = width;
	render_canvas.height = height;

	var text_canvas = this.text_canvas;
	text_canvas.width  = width;
	text_canvas.height = height

	this.camera.aspect = width / height;
	this.camera.updateProjectionMatrix();
	this.camControls.update();

	this.pathtracer.resize(width, height);
}

Snelly.prototype.resize = function()
{
	if (this.auto_resize)
	{
		let width = window.innerWidth;
		let height = window.innerHeight;
		this._resize(width, height);

		if (this.initialized)
			this.render();
	}
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
	// Check whether user is trying to click the Snelly home link, or user link
	var textCtx = this.textCtx;
	var x = event.pageX;
    var y = event.pageY;
	let linkWidth = this.textCtx.measureText('Snelly renderer vX.X.X').width;
	if (x>14 && x<14+linkWidth && y>15 && y<25) this.onSnellyLink = true;
	else this.onSnellyLink = false;
	if (this.sceneURL != '')
	{
		linkWidth = this.textCtx.measureText(this.sceneURL).width;
		if (x>14 && x<14+linkWidth && y>this.height-45 && y<this.height-35) this.onUserLink = true;
		else this.onUserLink = false;
	}

	this.camControls.update();
}

Snelly.prototype.onDocumentMouseDown = function(event)
{
	this.camControls.update();
}

Snelly.prototype.onDocumentMouseUp = function(event)
{
	this.camControls.update();
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
			this.camera.position.copy(this.initial_camera_position);
			this.camControls.target.copy(this.initial_camera_target);
			this.reset(true);
			break;

		case 82: // R key: reset scene 
			this.initScene();
			break;

		case 72: // H key: toggle hide/show dat gui
			this.guiVisible = !this.guiVisible;
			snelly.getGUI().toggleHide();
			break;
		
		case 79: // O key: output scene settings code to console
			let code = this.dumpScene();
			console.log(code);

		//case 87: console.log('w pressed'); break;
		//case 65: console.log('a pressed'); break;
		//case 83: console.log('s pressed'); break;
		//case 68: console.log('d pressed'); break;
	}
}

function camChanged()
{
	if (!snelly.rendering)
	{
		var no_recompile = true;
		snelly.reset(no_recompile);
	}
}





