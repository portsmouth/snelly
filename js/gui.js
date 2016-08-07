

var GUI = function() 
{
	// Create dat gui
	this.gui = new dat.GUI();
	var gui = this.gui;
	
	this.createLightTracerSettings();
	this.createSurfaceRendererSettings();
	this.createSceneSettings();
	this.createMaterialSettings();
	this.createEmissionSettings();
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
	this.lightTracerFolder.add(this.lightTracerSettings, 'exposure', -6.0, 6.0, 0.01);
	this.lightTracerFolder.add(this.lightTracerSettings, 'gamma', 0.0, 2.0, 0.01);
	this.lightTracerFolder.add(this.lightTracerSettings, 'maxPathLength', 1, 1024).onChange( function(value) { lightTracer.maxPathLength = Math.floor(value); lightTracer.reset(); } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'maxMarchSteps', 1, 1024).onChange( function(value) { lightTracer.maxMarchSteps = Math.floor(value); lightTracer.reset(); } );
	this.lightTracerFolder.add(this.lightTracerSettings, 'rayBufferSize', 16, 1024).onChange( function(value) { lightTracer.raySize = Math.floor(value); 
																											    lightTracer.initStates(); 
																												lightTracer.reset(); } );
	this.gui.remember(this.lightTracerSettings);
	this.lightTracerFolder.open();
}

GUI.prototype.createSurfaceRendererSettings = function()
{
	this.surfaceRendererFolder = this.gui.addFolder('Surface Renderer');
	this.surfaceRendererSettings = {};
	var surfaceRenderer = snelly.getSurfaceRenderer();

	this.surfaceRendererSettings.showSurface = surfaceRenderer.showSurface;
	this.surfaceRendererSettings.surfaceAlpha = surfaceRenderer.surfaceAlpha;
	this.surfaceRendererSettings.maxMarchSteps = surfaceRenderer.maxMarchSteps;
	
	this.surfaceRendererFolder.add(this.surfaceRendererSettings, 'showSurface').onChange( function(value)            { surfaceRenderer.showSurface = value;  } );
	this.surfaceRendererFolder.add(this.surfaceRendererSettings, 'surfaceAlpha', 0.0, 1.0).onChange( function(value) { surfaceRenderer.surfaceAlpha = value;  } );
	this.surfaceRendererFolder.add(this.surfaceRendererSettings, 'maxMarchSteps', 1, 1024).onChange( function(value) { surfaceRenderer.maxMarchSteps = Math.floor(value); surfaceRenderer.reset(); } );

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

	this.emissionFolder.add(this.emissionSettings, 'showLaserPointer').onChange( function(value) { laser.toggleVisibility(value); } );
	this.emissionFolder.add(laser, 'emissionRadius', 0.0, 10.0).onChange( function(value) { laser.setEmissionRadius(value);      lightTracer.reset(); } );
	this.emissionFolder.add(laser, 'emissionSpread', 0.0, 45.0).onChange( function(value) { laser.setEmissionSpreadAngle(value); lightTracer.reset(); } );
	this.emissionFolder.add(laser, 'emissionPower', 0.0, 10.0).onChange( function(value)  { laser.setEmissionPower(value);       lightTracer.reset(); } );
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




