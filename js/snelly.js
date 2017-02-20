


var Snelly = function()
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
	{
		this.stats = new Stats();
		this.stats.domElement.style.position = 'absolute';
		this.stats.domElement.style.left  = '0px'
		this.stats.domElement.style.bottom    = '0px'
		this.container.appendChild( this.stats.domElement );
	}

	// Text canvas
	var text_canvas = document.getElementById("text-canvas");
	this.textCtx = text_canvas.getContext("2d");
	this.onSnellyLink = false;

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

	// Setup laser pointer
	this.laser = new LaserPointer(this.glRenderer, this.glScene, this.camera, this.controls);
	this.laser.setPosition(new THREE.Vector3(-5.0, 0.0, 0.0));
	this.laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));

	// Instantiate scenes
	this.scenes = {}
	this.sceneObj = null;
	{
		this.addScene(new GemScene("Gem stone", ""));
		this.addScene(new ConvergingLensScene("Converging lens", ""));
		this.addScene(new DivergingLensScene("Diverging lens", ""));
		this.addScene(new LatticeScene("Lattice", ""));
		this.addScene(new SuperEllipsoidScene("Super-ellipsoid", ""));
		this.addScene(new FibreScene("Coiled fibre", ""));
		this.addScene(new TwistedTubeScene("Twisted tube", ""));
		this.addScene(new StackScene("Slab stack", ""));
		this.addScene(new TumblerScene("Tumbler", ""));
		this.addScene(new WineGlassScene("Wine glass", ""));
		this.addScene(new PoolScene("Pool", ""));
		this.addScene(new MengerScene("Menger sponge", ""));
		this.addScene(new KnotScene("Knot", ""));
		this.addScene(new MetaShapesScene("Metashapes", ""));
		this.addScene(new KIFSScene("KIFS fractal", ""));
		this.addScene(new WavesScene("Waves", ""));
		this.addScene(new PrismScene("Prism", ""));
		
		// ...
	}

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

		// Gases
		/*
		this.addMaterial( new Gas("Air", "", [0.0, 0.05792105, 238.0185, 0.00167917, 57.362]) );
		this.addMaterial( new Gas("Helium gas", "", [0.0, 0.01470091, 423.98]) );
		this.addMaterial( new Gas("Nitrogen gas", "", [6.497378e-5, 3.0738649e-2, 144.0]) );
		this.addMaterial( new Gas("Oxygen gas", "", [1.181494e-4, 9.708931e-3, 75.4]) );
		this.addMaterial( new Gas("Ammonia gas", "", [0.0, 0.032953, 90.392]) );
		this.addMaterial( new Gas("Argon gas", "", [0.0, 2.50141e-3, 91.012, 5.00283e-4, 87.892, 5.22343e-2, 214.02]) );
		this.addMaterial( new Gas("Neon gas", "", [0.0, 0.00128145, 184.661, 0.0220486, 376.840]) );
		this.addMaterial( new Gas("Krypton gas", "", [0.0, 0.00253637, 65.4742, 0.00273649, 73.698, 0.0620802, 181.08]) );
		this.addMaterial( new Gas("Xenon gas", "", [0.0, 0.00322869, 46.301, 0.00355393, 50.578, 0.0606764, 112.74]) );
		*/

		// @todo: Metals. Need a much better approximation than 'LinearMetal' for plausible metals 
		this.addMetal( new LinearMetal("Aluminium", "",  0.46555, 4.7121, 1.6620, 8.0439) );
		this.addMetal( new LinearMetal("Gold", "",       1.5275, 1.8394, 0.16918, 3.8816) );
	}

	// Instantiate light tracer
	this.lightTracer = new LightTracer();

	// Instantiate distance field surface renderer
	this.surfaceRenderer = new SurfaceRenderer();
	
	// Do initial resize:
	this.resize();

	// Load the initial scene and material
	this.gui = null;

	this.loadDielectric("Glass (LASF35)");
	this.loadMetal("Gold");

	this.loadScene("Gem stone");

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

			case 70: // F key: focus on emitter
				snelly.controls.object.zoom = snelly.controls.zoom0;
				snelly.controls.target.copy(snelly.laser.getPoint());
				snelly.controls.update();
				break;

			case 82: // R key: reset scene 
				// @todo
				break;

			case 72: // H key: hide dat gui
				snelly.getGUI().visible = !(snelly.getGUI().visible);
				if  (snelly.getGUI().visible) { snelly.stats.domElement.style.visibility = "visible"; }
				else                          { snelly.stats.domElement.style.visibility = "hidden"; }
				break;

			case 67: // C key: dev tool to dump cam and laser details, for setting scene defaults
				var lp = snelly.laser.getPosition();
				var ld = snelly.laser.getDirection();
				var er = snelly.laser.getEmissionRadius();
				var ex = snelly.laser.getEmissionSpreadAngle();
				var t = snelly.controls.target;
				var c = snelly.camera.position;

				console.log(`laser.setPosition(new THREE.Vector3(${lp.x.toPrecision(6)}, ${lp.y.toPrecision(6)}, ${lp.z.toPrecision(6)}));`);
				if (snelly.laser.targetSet)
				{
					var tg = snelly.laser.target;
					console.log(`laser.setTarget(new THREE.Vector3(${tg.x.toPrecision(6)}, ${tg.y.toPrecision(6)}, ${tg.z.toPrecision(6)}));`);
				}
				else
				{
					console.log(`laser.setDirection(new THREE.Vector3(${ld.x.toPrecision(6)}, ${ld.y.toPrecision(6)}, ${ld.z.toPrecision(6)}));`);
				}
				console.log(`laser.setEmissionRadius(${er.toPrecision(6)});`);
				console.log(`laser.setEmissionSpreadAngle(${ex.toPrecision(6)});`);
				console.log(`controls.target.set(${t.x.toPrecision(6)}, ${t.y.toPrecision(6)}, ${t.z.toPrecision(6)});`);
				console.log(`camera.position.set(${c.x.toPrecision(6)}, ${c.y.toPrecision(6)}, ${c.z.toPrecision(6)});`);
				break;
			
		}
	}, false);

	this.glRenderer.domElement.addEventListener( 'mousemove', this, false );
	this.glRenderer.domElement.addEventListener( 'mousedown', this, false );
	this.glRenderer.domElement.addEventListener( 'mouseup',   this, false );
	this.glRenderer.domElement.addEventListener( 'contextmenu',   this, false );
	this.glRenderer.domElement.addEventListener( 'click', this, false );

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
	this.laser.unsetTarget();

	this.sceneObj = this.scenes[sceneName];
	this.sceneObj.init(this.controls, this.camera, this.laser);
	this.laser.buildEmitterGeo();

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

Snelly.prototype.getLoadedScene = function()
{
	return this.sceneObj;
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
Snelly.prototype.reset = function()
{	
	this.surfaceRenderer.reset();
	this.lightTracer.reset();

	if (!this.initialized) return;

	this.gui.sync();

	this.render();
}

// Render all 
Snelly.prototype.render = function()
{
	if (!this.initialized) return;
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

	// Update HUD text canvas	
	this.textCtx.textAlign = "left";   	// This determines the alignment of text, e.g. left, center, right
	this.textCtx.textBaseline = "middle";	// This determines the baseline of the text, e.g. top, middle, bottom
	this.textCtx.font = '12px monospace';	// This determines the size of the text and the font family used
	
	this.textCtx.clearRect(0, 0, this.textCtx.canvas.width, this.textCtx.canvas.height);
	this.textCtx.globalAlpha = 0.95;

	if (snelly.getGUI().visible)
	{
	  	var lsStats = this.lightTracer.getStats();

	  	if (this.onSnellyLink)
	  	{
	  		this.textCtx.fillStyle = "#ff5500";
	  	}
	  	else
	  	{
	  		this.textCtx.fillStyle = "#ffff00";
	  	}
	  	
	  	this.textCtx.fillText('Snelly renderer', 14, 20);

	  	this.textCtx.fillStyle = "#aaaaff";
	  	this.textCtx.fillText('ray count:    ' + (lsStats.rayCount/1.0e6).toPrecision(3) + 'M', 14, 35);
	  	this.textCtx.fillText('waves traced: ' + lsStats.wavesTraced,    14, 50);
	}

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

	var text_canvas = document.getElementById("text-canvas");
	text_canvas.width  = width;
	text_canvas.height = height

	this.camera.aspect = width / height;
	this.camera.updateProjectionMatrix();

	this.lightTracer.resize(width, height);
	this.surfaceRenderer.resize(width, height);
	this.glRenderer.setSize(width, height);
	this.laser.resize(width, height);

	if (this.initialized)
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
		case 'click':       this.onClick(event);  break;
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
	// Check whether user is trying to click the Snelly home link
	var textCtx = this.textCtx;
	var x = event.layerX;
    var y = event.layerY;
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
	this.controls.update();
	event.preventDefault();
	if (event.altKey) return; // don't pick if alt-right-clicking (panning)

	var xPick =   (( event.clientX - this.glRenderer.domElement.offsetLeft ) / this.glRenderer.domElement.width)*2 - 1;
	var yPick = - (( event.clientY - this.glRenderer.domElement.offsetTop ) / this.glRenderer.domElement.height)*2 + 1;

	var pickedPoint = this.surfaceRenderer.pick(xPick, yPick);
	if (pickedPoint == null)
	{
		this.laser.unsetTarget();
		return;
	} 

	this.laser.setTarget(pickedPoint);	
	this.reset();
}

function camChanged()
{
	snelly.reset();
	snelly.render();
}





