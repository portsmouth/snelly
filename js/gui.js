

var GUI = function() 
{
	// Create dat gui
	this.gui = new dat.GUI();
	this.gui.domElement.id = 'gui';
	var gui = this.gui;

	this.visible = true;
	
	this.createSceneSettings();
	this.createMaterialSettings();
	this.createEmissionSettings();
	this.createLightTracerSettings();
	this.createPathtracerSettings();
}

function updateDisplay(gui) {
    for (var i in gui.__controllers) {
        gui.__controllers[i].updateDisplay();
    }
    for (var f in gui.__folders) {
        updateDisplay(gui.__folders[f]);
    }
}

GUI.prototype.sync = function()
{
	var laser = snelly.getLaser();
	this.emissionSettings.eulerAngles.x = laser.eulerAngles.x * 180.0/Math.PI;
	this.emissionSettings.eulerAngles.y = laser.eulerAngles.y * 180.0/Math.PI;
	this.emissionSettings.eulerAngles.z = laser.eulerAngles.z * 180.0/Math.PI;
	this.emissionSettings.emissionRadius = laser.getEmissionRadius();
	this.emissionSettings.emissionSpread = laser.getEmissionSpreadAngle();

	updateDisplay(this.gui);
}

GUI.prototype.createLightTracerSettings = function()
{
	this.lightTracerFolder = this.gui.addFolder('Light Tracer');
	this.lightTracerSettings = {};
	var lightTracer = snelly.getLightTracer();

	this.lightTracerSettings.enable = lightTracer.enabled;
	this.lightTracerSettings.exposure = 0.0;
	this.lightTracerSettings.gamma = 2.2;
	this.lightTracerSettings.maxPathLength = lightTracer.maxPathLength;
	this.lightTracerSettings.maxMarchSteps = lightTracer.maxMarchSteps;
	this.lightTracerSettings.rayBufferSize = lightTracer.raySize;
	
	this.lightTracerFolder.add(this.lightTracerSettings, 'enable').onChange( function(value) { lightTracer.enabled = value;  } );
	this.lightTracerFolder.add(lightTracer, 'showExternal');
	this.lightTracerFolder.add(lightTracer, 'showInternal');
	this.lightTracerFolder.add(this.lightTracerSettings, 'exposure', -10.0, 10.0, 0.01);
	this.lightTracerFolder.add(this.lightTracerSettings, 'gamma', 0.0, 4.0, 0.01);
	this.lightTracerFolder.add(this.lightTracerSettings, 'maxPathLength', 4, 1024).onChange( function(value) { lightTracer.maxPathLength = Math.floor(value); lightTracer.reset(); } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'maxMarchSteps', 32, 1024).onChange( function(value) { lightTracer.maxMarchSteps = Math.floor(value); lightTracer.reset(); } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'rayBufferSize', 64, 1024).onChange( function(value) { lightTracer.raySize = Math.floor(value); 
																											    lightTracer.initStates(); 
																												lightTracer.reset(); } );
	this.gui.remember(this.lightTracerSettings);
	this.lightTracerFolder.close();
}


function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}


GUI.prototype.createPathtracerSettings = function()
{
	this.pathtracerFolder = this.gui.addFolder('Pathtracer');
	this.pathtracerSettings = {};
	var pathtracer = snelly.getPathtracer();

	var renderModes = ['normals', 'blinn'];
	
	this.pathtracerFolder.add(pathtracer, 'enable');
	//this.pathtracerFolder.add(pathtracer, 'depthTest');
	this.pathtracerFolder.add(pathtracer, 'showBounds');
	//this.pathtracerFolder.add(pathtracer, 'renderMode', renderModes).onChange( function(renderMode) { pathtracer.reset(); });
	this.pathtracerFolder.add(pathtracer, 'exposure', 0.0, 50.0);
	this.pathtracerFolder.add(pathtracer, 'maxBounces', 1, 10).onChange( function(value) { pathtracer.maxBounces = Math.floor(value); pathtracer.reset(); } );
	this.pathtracerFolder.add(pathtracer, 'maxMarchSteps', 1, 1024).onChange( function(value) { pathtracer.maxMarchSteps = Math.floor(value); pathtracer.reset(); } );

	this.gui.remember(this.pathtracerSettings);
	this.pathtracerFolder.open();
}


GUI.prototype.createEmissionSettings = function()
{
	this.emissionFolder = this.gui.addFolder('Emission');
	var lightTracer = snelly.getLightTracer();
	var pathtracer = snelly.getPathtracer();
	var laser = snelly.getLaser();
	this.emissionSettings = {};
	this.emissionSettings.showLaserPointer = false;
	this.emissionSettings.spectrum = 'monochromatic';
	this.emissionSettings.emissionRadius = 0.01;

	this.emissionFolder.add(this.emissionSettings, 'showLaserPointer').onChange( function(value) { laser.toggleVisibility(value); } );
	
	this.emissionRadiusControl = this.emissionFolder.add(this.emissionSettings, 'emissionRadius', 0.0, 3.0);
	this.emissionRadiusControl.onChange( function(value) 
	{ 
		laser.setEmissionRadius(value);      
		lightTracer.reset(); 
		pathtracer.reset();
	} );
	this.emissionFolder.add(laser, 'emissionSpread', 0.0, 90.0).onChange( function(value) 
	{ 
		laser.setEmissionSpreadAngle(value);
		lightTracer.reset(); 
		pathtracer.reset();
	} );
	this.emissionFolder.add(laser, 'emissionPower', 0.0, 10.0).onChange( function(value) 
	{ 
		laser.setEmissionPower(value);
		lightTracer.reset(); 
		pathtracer.reset();
	} );

	this.emissionFolder.add(laser, 'skyPower', 0.0, 1.0).onChange( function(value) 
	{ 
		laser.setSkyPower(value);
		lightTracer.reset(); 
		pathtracer.reset();
	} );

	this.gui.remember(laser);

	// Spectrum selection
	var GUI = this;
	var spectrumObj = snelly.getLoadedSpectrum();
	var spectrumName = spectrumObj.getName();
	var spectra = snelly.getSpectra();
	var spectrumNames = Object.keys(spectra);

	this.emissionSettings["spectrum selection"] = spectrumName;
	this.emissionFolder.add(this.emissionSettings, 'spectrum selection', spectrumNames).onChange( function(spectrumName) {

						// remove gui for current spectrum
						var spectrumObj = snelly.getLoadedSpectrum();
						spectrumObj.eraseGui(GUI.emissionFolder);

						// load new spectrum
				 		snelly.loadSpectrum(spectrumName);
	
				 		// init gui for new spectrum
				 		spectrumObj = snelly.getLoadedSpectrum();
				 		spectrumObj.initGui(GUI.emissionFolder);
				 		
				 	} );

	spectrumObj.initGui(this.emissionFolder);

	// Laser transform sliders
	this.emissionTranslationFolder = this.emissionFolder.addFolder('Translation');
	var xT = this.emissionTranslationFolder.add(laser.group.position, 'x');
	var yT = this.emissionTranslationFolder.add(laser.group.position, 'y');
	var zT = this.emissionTranslationFolder.add(laser.group.position, 'z');
	xT.onChange( function(value)  { snelly.reset(); } );
	yT.onChange( function(value)  { snelly.reset(); } );
	zT.onChange( function(value)  { snelly.reset(); } );

	this.emissionSettings.eulerAngles = new THREE.Vector3();
	this.emissionSettings.eulerAngles.x = laser.eulerAngles.x * 180.0/Math.PI;
	this.emissionSettings.eulerAngles.y = laser.eulerAngles.y * 180.0/Math.PI;
	this.emissionSettings.eulerAngles.z = laser.eulerAngles.z * 180.0/Math.PI;

	this.emissionRotationFolder = this.emissionFolder.addFolder('Rotation');
	var xR = this.emissionRotationFolder.add(this.emissionSettings.eulerAngles, 'x', -180.0, 180.0);
	xR.onChange( function(value) 
	{ 
		var euler = laser.getEuler();
		euler.x = value * Math.PI/180.0;
		laser.setEuler(euler);
		snelly.reset();
	} );
	var yR = this.emissionRotationFolder.add(this.emissionSettings.eulerAngles, 'y', -180.0, 180.0);
	yR.onChange( function(value) 
	{ 
		var euler = laser.getEuler();
		euler.y = value * Math.PI/180.0;
		laser.setEuler(euler);
		snelly.reset();
	} );
	var zR = this.emissionRotationFolder.add(this.emissionSettings.eulerAngles, 'z', -180.0, 180.0);
	zR.onChange( function(value) 
	{ 
		var euler = laser.getEuler();
		euler.z = value * Math.PI/180.0;
		laser.setEuler(euler);
		snelly.reset();
	} );

	this.gui.remember(this.emissionSettings);
	this.emissionFolder.open();
	return this.emissionFolder;
}


GUI.prototype.createSceneSettings = function()
{
	this.sceneFolder = this.gui.addFolder('Scene');
	var sceneObj = snelly.getLoadedScene();
	var sceneName = sceneObj.getName();
	scenes = snelly.getScenes();
	var sceneNames = Object.keys(scenes);

	// Scene selection menu
	this.sceneSettings = {};
	this.sceneSettings["scene selection"] = sceneName;
	var GUI = this;

	this.sceneFolder.add(this.sceneSettings, 'scene selection', sceneNames).onChange( function(sceneName) {

						// remove gui for current scene
						var sceneObj = snelly.getLoadedScene();
						sceneObj.eraseGui(GUI.sceneFolder);

						// load new scene
				 		snelly.loadScene(sceneName);

				 		// init gui for new scene
				 		sceneObj = snelly.getLoadedScene();
				 		sceneObj.initGui(GUI.sceneFolder);
				 		
				 	} );

	sceneObj.initGui(this.sceneFolder);
	this.sceneFolder.close();

	this.gui.remember(this.sceneSettings);
	return this.sceneFolder;
}


GUI.prototype.createMaterialSettings = function()
{
	var GUI = this;
	this.materialFolder = this.gui.addFolder('Material');

	// Dielectric settings
	this.dielectricFolder = this.materialFolder.addFolder('Dielectric');
	var dielectricObj = snelly.getLoadedDielectric();
	var dielectricName = dielectricObj.getName();
	dielectrics = snelly.getDielectrics();
	var dielectricNames = Object.keys(dielectrics);
	this.dielMaterialSettings = {};
	this.dielMaterialSettings["dielectric material"] = dielectricName;
	this.dielectricFolder.add(this.dielMaterialSettings, 'dielectric material', dielectricNames).onChange( function(materialName) {
					
						// remove gui for current material
						var materialObj = snelly.getLoadedDielectric();
						materialObj.eraseGui(GUI.dielectricFolder);
						
						// load new material
				 		snelly.loadDielectric(materialName);
				 		
				 		// init gui for new material
				 		materialObj = snelly.getLoadedDielectric();
				 		materialObj.initGui(GUI.dielectricFolder);
				 	} );
	dielectricObj.initGui(this.dielectricFolder);

	// Metal settings
	this.metalFolder = this.materialFolder.addFolder('Metal');
	var metalObj = snelly.getLoadedMetal();
	var metalName = metalObj.getName();
	metals = snelly.getMetals();
	var metalNames = Object.keys(metals);
	this.metalMaterialSettings = {};
	this.metalMaterialSettings["metal material"] = metalName;
	this.metalFolder.add(this.metalMaterialSettings, 'metal material', metalNames).onChange( function(materialName) {
						
						// remove gui for current material
						var materialObj = snelly.getLoadedMetal();
						materialObj.eraseGui(GUI.metalFolder);
						
						// load new material
				 		snelly.loadMetal(materialName);
				 		
				 		// init gui for new material
				 		materialObj = snelly.getLoadedMetal();
				 		materialObj.initGui(GUI.metalFolder);
			 	} );
	metalObj.initGui(this.metalFolder);
	
	this.materialFolder.open();
	this.gui.remember(this.materialSettings);
	return this.materialFolder;
}




