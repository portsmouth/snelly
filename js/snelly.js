
var Snelly = function()
{
	this.initialized = false; 

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;

	window.addEventListener( 'resize', this, false );

	this.container = document.getElementById('container');
	{
		this.stats = new Stats();
		this.stats.domElement.style.position = 'absolute';
		this.stats.domElement.style.top = '0px';
		this.container.appendChild( this.stats.domElement );
	}

	// Setup THREE.js GL viewport renderer and camera
	var ui_canvas = document.getElementById('ui-canvas');
	ui_canvas.style.top = 0;
	ui_canvas.style.position = 'fixed' 

	var VIEW_ANGLE = 45;
	var ASPECT = this.width / this.height ;
	var NEAR = 0.05;
	var FAR = 1000;
	this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);

	this.glRenderer = new THREE.WebGLRenderer( { canvas: ui_canvas,
											     alpha: true,
											     antialias: true } );
	this.glRenderer.setClearColor( 0x000000, 0 ); // the default
	this.glRenderer.setSize(this.width, this.height);
	this.glScene = new THREE.Scene();
	this.glScene.add(this.camera);

	var pointLight = new THREE.PointLight(0xa0a0a0);
	pointLight.position.x = 10;
	pointLight.position.y = 50;
	pointLight.position.z = 130;
	this.glScene.add(pointLight);

	var light = new THREE.AmbientLight( 0x808080 ); // soft white light
	this.glScene.add( light );


	// Create user control system for camera
	this.controls = new THREE.OrbitControls(this.camera, this.glRenderer.domElement);
	this.controls.zoomSpeed = 2.0;
	this.controls.addEventListener( 'change', camChanged );

	// Setup Laser pointer
	this.laser = new LaserPointer(this.glRenderer, this.glScene, this.camera, this.controls);
	this.laser.setPosition(new THREE.Vector3(-5.0, 0.0, 0.0));
	this.laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));

	// Setup keypress and mouse events
	snelly = this;
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
			case 70: // F key: focus on emitter
				snelly.controls.object.zoom = snelly.controls.zoom0;
				snelly.controls.target.copy(snelly.laser.getPoint());
				snelly.controls.update();
				break; 
		}
	}, false);
	
	this.glRenderer.domElement.addEventListener( 'mousemove', this, false );
	this.glRenderer.domElement.addEventListener( 'mousedown', this, false );
	this.glRenderer.domElement.addEventListener( 'mouseup',   this, false );
	this.glRenderer.domElement.addEventListener( 'contextmenu',   this, false );

	// Instantiate scenes
	this.scenes = {}
	this.sceneObj = null;
	{
		this.addScene(new SphereScene("sphere", "Simple sphere"));
		this.addScene(new FibreScene("fibre", "Simple optical fibre"));
		this.addScene(new BoxScene("box", "Simple box"));
		this.addScene(new OceanScene("ocean", "Ocean"));
		this.addScene(new MengerScene("menger", "Menger Sponge"));
		this.addScene(new EllipsoidScene("ellipsoid", "Ellipsoid"));

		// ...
	}

	// Instantiate materials
	this.materials = {}
	this.materialObj = null;
	{
		// Dielectrics
		this.addMaterial( new ConstantDielectric("constant", "Constant IOR dielectric", 1.5) ); 
		this.addMaterial( new SellmeierDielectric("diamond", "Diamond",       0.0, 0.3306,     0.175,         4.3356,      0.1060,       0.0,        0.0       ) );
		
		var glass = new SellmeierDielectric("glass_bk7", "Glass (BK7)", 0.0, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653);
		this.addMaterial(glass);

		// Metals
		this.addMaterial( new LinearMetal("aluminium", "Aluminium",  0.46555, 4.7121, 1.6620, 8.0439) );
		this.addMaterial( new LinearMetal("gold", "Gold",            1.5275, 1.8394, 0.16918, 3.8816) );
	
		// ...
	}

	// Instantiate light tracer
	this.lightTracer = new LightTracer();

	// Instantiate distance field surface renderer
	this.surfaceRenderer = new SurfaceRenderer();

	// Do initial resize:
	this.resize();

	// Load the initial scene and material
	this.loadScene("sphere");
	this.loadMaterial("glass_bk7");

	// Create dat gui
	this.gui = new GUI();

	this.initialized = true; 
}

Snelly.prototype.getLightTracer = function()
{
	return this.lightTracer;
}

Snelly.prototype.getSurfaceRenderer = function()
{
	return this.surfaceRenderer;
}

Snelly.prototype.getGUI= function()
{
	return this.gui;
}

Snelly.prototype.getLaser = function()
{
	return this.laser;
}

Snelly.prototype.getCamera = function()
{
	return this.camera;
}


//
// Scene management
//
Snelly.prototype.addScene = function(sceneObj)
{
	this.scenes[sceneObj.getName()] = sceneObj;
}

Snelly.prototype.getScenes = function()
{
	return this.scenes;
}

Snelly.prototype.loadScene = function(sceneName)
{
	this.sceneObj = this.scenes[sceneName];
	this.sceneObj.setCam(this.controls, this.camera);
	this.sceneObj.setLaser(this.laser);

	// Camera frustum update
	this.camera.near = 1.0e-2*this.sceneObj.getScale();
	this.camera.far  = 1.0e2*this.sceneObj.getScale();

	this.reset();
}

Snelly.prototype.getLoadedScene = function()
{
	return this.sceneObj;
}


//
// Material management
//
Snelly.prototype.addMaterial = function(materialObj)
{
	this.materials[materialObj.getName()] = materialObj;
}

Snelly.prototype.getMaterials = function()
{
	return this.materials;
}

Snelly.prototype.loadMaterial = function(materialName)
{
	this.materialObj = this.materials[materialName];
	this.reset();
}

Snelly.prototype.getLoadedMaterial = function()
{
	return this.materialObj;
}


var flag = 0;

// Renderer reset on camera update
Snelly.prototype.reset = function()
{	
	this.lightTracer.reset();
	this.surfaceRenderer.reset();
	if (this.initialized)
		this.render();
}

// Render all 
Snelly.prototype.render = function()
{
	if (this.sceneObj == null) return;
	
	var gl = GLU.gl;
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.enable(gl.DEPTH_TEST);

	gl.viewport(0, 0, this.width, this.height);

	// Render laser pointer
	this.laser.render();

	// Render light beams
	this.lightTracer.render();

	// Render distance field surface
	this.surfaceRenderer.render();

	// Update stats
	this.stats.update();
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

	var ui_canvas = document.getElementById('ui-canvas');
	ui_canvas.width  = width;
	ui_canvas.height = height;

	this.camera.aspect = width / height;
	this.camera.updateProjectionMatrix();

	this.lightTracer.resize(width, height);
	this.surfaceRenderer.resize(width, height);
	this.glRenderer.setSize(width, height);

	this.render();
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
	}
}

Snelly.prototype.onDocumentMouseMove = function(event)
{
	this.controls.update();
	event.preventDefault();
	if (this.laser.onMouseMove(event)) this.reset();
}

Snelly.prototype.onDocumentMouseDown = function(event)
{
	this.controls.update();
	event.preventDefault();
	this.laser.onMouseDown(event);
}

Snelly.prototype.onDocumentMouseUp = function(event)
{
	this.controls.update();
	event.preventDefault();
	this.laser.onMouseUp(event);
}

Snelly.prototype.onDocumentRightClick = function(event)
{
	// @todo: instead, do pick on left click,
	//        rotate on alt-left click, save right click for pan

	this.controls.update();
	event.preventDefault();

	var xPick =   (( event.clientX - this.glRenderer.domElement.offsetLeft ) / this.glRenderer.domElement.width)*2 - 1;
	var yPick = - (( event.clientY - this.glRenderer.domElement.offsetTop ) / this.glRenderer.domElement.height)*2 + 1;

	var pickedPoint = this.surfaceRenderer.pick(xPick, yPick);
	if (pickedPoint == null) return;

	this.laser.setTarget(pickedPoint);	
	this.reset();
}

function camChanged()
{
	snelly.reset();
	snelly.render();
}





