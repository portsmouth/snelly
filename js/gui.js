

var GUI = function() 
{
	// Create dat gui
	this.gui = new dat.GUI();
	var gui = this.gui;
	
	this.createSceneSettings();
	this.createMaterialSettings();
	this.createEmissionSettings();
	this.createLightTracerSettings();
	this.createSurfaceRendererSettings();
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

	updateDisplay(this.gui);
}

GUI.prototype.createLightTracerSettings = function()
{
	this.lightTracerFolder = this.gui.addFolder('Light Tracer');
	this.lightTracerSettings = {};
	var lightTracer = snelly.getLightTracer();

	this.lightTracerSettings.enable = lightTracer.enabled;
	this.lightTracerSettings.exposure = 1.0;
	this.lightTracerSettings.gamma = 2.2;
	this.lightTracerSettings.maxPathLength = lightTracer.maxPathLength;
	this.lightTracerSettings.maxMarchSteps = lightTracer.maxMarchSteps;
	this.lightTracerSettings.rayBufferSize = lightTracer.raySize;
	
	this.lightTracerFolder.add(this.lightTracerSettings, 'enable').onChange( function(value) { lightTracer.enabled = value;  } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'exposure', -2.0, 4.0, 0.01);
	this.lightTracerFolder.add(this.lightTracerSettings, 'gamma', 0.0, 4.0, 0.01);
	this.lightTracerFolder.add(this.lightTracerSettings, 'maxPathLength', 4, 1024).onChange( function(value) { lightTracer.maxPathLength = Math.floor(value); lightTracer.reset(); } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'maxMarchSteps', 32, 1024).onChange( function(value) { lightTracer.maxMarchSteps = Math.floor(value); lightTracer.reset(); } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'rayBufferSize', 64, 1024).onChange( function(value) { lightTracer.raySize = Math.floor(value); 
																											    lightTracer.initStates(); 
																												lightTracer.reset(); } );
	this.gui.remember(this.lightTracerSettings);
	this.lightTracerFolder.open();
}


function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}


GUI.prototype.createSurfaceRendererSettings = function()
{
	this.surfaceRendererFolder = this.gui.addFolder('Surface Renderer');
	this.surfaceRendererSettings = {};
	var surfaceRenderer = snelly.getSurfaceRenderer();

	var renderModes = ['normals', 'blinn'];
	
	this.surfaceRendererFolder.add(surfaceRenderer, 'enable');
	this.surfaceRendererFolder.add(surfaceRenderer, 'depthTest');
	this.surfaceRendererFolder.add(surfaceRenderer, 'showBounds');

	this.surfaceRendererFolder.add(surfaceRenderer, 'renderMode', renderModes).onChange( function(renderMode) { surfaceRenderer.reset(); });
	this.surfaceRendererFolder.add(surfaceRenderer, 'surfaceAlpha', 0.0, 1.0);
	this.surfaceRendererFolder.add(surfaceRenderer, 'maxMarchSteps', 1, 1024).onChange( function(value) { surfaceRenderer.maxMarchSteps = Math.floor(value); surfaceRenderer.reset(); } );

	this.surfaceRendererSettings.diffuseCol1 = [surfaceRenderer.kd1[0]*255.0, surfaceRenderer.kd1[1]*255.0, surfaceRenderer.kd1[2]*255.0];
	this.surfaceRendererSettings.diffuseCol2 = [surfaceRenderer.kd2[0]*255.0, surfaceRenderer.kd2[1]*255.0, surfaceRenderer.kd2[2]*255.0];

	this.surfaceRendererFolder.addColor(this.surfaceRendererSettings, 'diffuseCol1').onChange( function(value) 
	{
		if (typeof value==='string' || value instanceof String)
		{
			var color = hexToRgb(value);
			surfaceRenderer.kd1[0] = color.r / 255.0;
			surfaceRenderer.kd1[1] = color.g / 255.0;
			surfaceRenderer.kd1[2] = color.b / 255.0;
		}
		else
		{
			surfaceRenderer.kd1[0] = value[0] / 255.0;
			surfaceRenderer.kd1[1] = value[1] / 255.0;
			surfaceRenderer.kd1[2] = value[2] / 255.0;
		}
		surfaceRenderer.reset(); 
	});

	this.surfaceRendererFolder.addColor(this.surfaceRendererSettings, 'diffuseCol2').onChange( function(value) 
	{ 
		if (typeof value==='string' || value instanceof String)
		{
			var color = hexToRgb(value);
			surfaceRenderer.kd2[0] = color.r / 255.0;
			surfaceRenderer.kd2[1] = color.g / 255.0;
			surfaceRenderer.kd2[2] = color.b / 255.0;
		}
		else
		{
			surfaceRenderer.kd2[0] = value[0] / 255.0;
			surfaceRenderer.kd2[1] = value[1] / 255.0;
			surfaceRenderer.kd2[2] = value[2] / 255.0;
		}
		surfaceRenderer.reset(); 
	});

	this.surfaceRendererFolder.add(surfaceRenderer, 'specPower', 1.0, 100.0).onChange( function(renderMode) { surfaceRenderer.reset(); });

	this.gui.remember(this.surfaceRendererSettings);
	this.surfaceRendererFolder.open();
}


GUI.prototype.createEmissionSettings = function()
{
	this.emissionFolder = this.gui.addFolder('Emission');
	var lightTracer = snelly.getLightTracer();
	var laser = snelly.getLaser();
	this.emissionSettings = {};
	this.emissionSettings.showLaserPointer = true;
	this.emissionSettings.spectrum = 'monochromatic';
	this.emissionSettings.emissionRadius = 0.01; // (scene-scale relative)

	this.emissionFolder.add(this.emissionSettings, 'showLaserPointer').onChange( function(value) { laser.toggleVisibility(value); } );
	this.emissionFolder.add(this.emissionSettings, 'emissionRadius', 0.0, 3.0).onChange( function(value) 
	{ 
		var sceneObj = snelly.getLoadedScene();
		var sceneScale = sceneObj.getScale();
		laser.setEmissionRadius(sceneScale*value);      
		lightTracer.reset(); 
	} );
	this.emissionFolder.add(laser, 'emissionSpread', 0.0, 45.0).onChange( function(value) { laser.setEmissionSpreadAngle(value); lightTracer.reset(); } );
	this.gui.remember(laser);

	// Spectrum selection
	var GUI = this;
	var spectrumObj = lightTracer.getLoadedSpectrum();
	var spectrumName = spectrumObj.getName();
	var spectra = lightTracer.getSpectra();
	var spectrumNames = Object.keys(spectra);

	this.emissionSettings["spectrum selection"] = spectrumName;
	this.emissionFolder.add(this.emissionSettings, 'spectrum selection', spectrumNames).onChange( function(spectrumName) {

						// remove gui for current spectrum
						var spectrumObj = lightTracer.getLoadedSpectrum();
						spectrumObj.eraseGui(GUI.emissionFolder);

						// load new scene
				 		lightTracer.loadSpectrum(spectrumName);

				 		// init gui for new scene
				 		spectrumObj = lightTracer.getLoadedSpectrum();
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
	this.sceneFolder.open();

	this.gui.remember(this.sceneSettings);
	return this.sceneFolder;
}


GUI.prototype.createMaterialSettings = function()
{
	this.materialFolder = this.gui.addFolder('Material');
	var materialObj = snelly.getLoadedMaterial();
	var materialName = materialObj.getName();
	materials = snelly.getMaterials();
	var materialNames = Object.keys(materials);

	// Material selection menu
	this.materialSettings = {};
	this.materialSettings["material selection"] = materialName;
	var GUI = this;

	this.materialFolder.add(this.materialSettings, 'material selection', materialNames).onChange( function(materialName) {

						// remove gui for current material
						var materialObj = snelly.getLoadedMaterial();
						materialObj.eraseGui(GUI.materialFolder);

						// load new material
				 		snelly.loadMaterial(materialName);

				 		// init gui for new material
				 		materialObj = snelly.getLoadedMaterial();
				 		materialObj.initGui(GUI.materialFolder);
				 		
				 	} );
	
	materialObj.initGui(this.materialFolder);
	this.materialFolder.open();

	this.gui.remember(this.materialSettings);
	return this.materialFolder;
}




