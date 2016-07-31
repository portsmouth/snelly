

var GUI = function(renderer) 
{
	this.renderer = renderer;
	
	// Create dat gui
	this.gui = new dat.GUI();
	var gui = this.gui;
	
	this.rs = this.rendererSettings();
	this.ss = this.sceneSettings();
	this.ss = this.materialSettings();
	this.es = this.emissionSettings();

	this.registerCamera();
	this.registerLaser();
}

GUI.prototype.rendererSettings = function()
{
	this.rendererFolder = this.gui.addFolder('Renderer');
	var renderer = this.renderer;
	renderer.renderSettings = {};
	settings = renderer.renderSettings;

	settings.exposure = 10.0;
	settings.showSurface = false;
	settings.surfaceAlpha = 0.5;
	settings.maxPathLength = renderer.maxPathLength;
	settings.maxMarchSteps =renderer.maxMarchSteps;
	settings.photonsPerFrame = renderer.raySize;
	
	this.rendererFolder.add(settings, 'exposure', 0.0, 10.0, 0.01);
	this.rendererFolder.add(settings, 'maxPathLength', 1, 1024).onChange( function(value) { renderer.maxPathLength = Math.floor(value); renderer.reset(); } );
	this.rendererFolder.add(settings, 'maxMarchSteps', 1, 1024).onChange( function(value) { renderer.maxMarchSteps = Math.floor(value); renderer.reset(); } );
	this.rendererFolder.add(settings, 'photonsPerFrame', 1, 1024).onChange( function(value) { renderer.raySize = Math.floor(value); renderer.initRayStates(); renderer.reset(); } );
	this.rendererFolder.add(settings, 'showSurface');
	this.rendererFolder.add(settings, 'surfaceAlpha');
	this.gui.remember(settings);

	this.rendererFolder.open();
	return this.rendererFolder;
}

GUI.prototype.registerCamera = function()
{
	/*
	rs = this.rs;
	var m = 100.0;
	this.cpfolder = rs.addFolder('Camera Position');
	this.cpfolder.add(this.renderer.camera.position, 'x', -m, m).listen();
	this.cpfolder.add(this.renderer.camera.position, 'y', -m, m).listen();
	this.cpfolder.add(this.renderer.camera.position, 'z', -m, m).listen();
	this.crfolder = rs.addFolder('Camera Target');
	this.crfolder.add(this.renderer.controls.target, 'x', -m, m).listen();
	this.crfolder.add(this.renderer.controls.target, 'y', -m, m).listen();
	this.crfolder.add(this.renderer.controls.target, 'z', -m, m).listen();
	this.gui.remember(this.renderer.camera.position);
	this.gui.remember(this.renderer.camera.rotation);
	this.cpfolder.close();
	this.crfolder.close();
	*/
}

GUI.prototype.registerLaser = function()
{
	/*
	es = this.es;
	var m = 7.0;
	this.lpfolder = es.addFolder('Laser Position');
	this.lpfolder.add(this.renderer.laser.group.position, 'x', -m, m).listen();
	this.lpfolder.add(this.renderer.laser.group.position, 'y', -m, m).listen();
	this.lpfolder.add(this.renderer.laser.group.position, 'z', -m, m).listen();
	this.lrfolder = es.addFolder('Laser Rotation');
	this.lrfolder.add(this.renderer.laser.group.rotation, 'x', -m, m).listen();
	this.lrfolder.add(this.renderer.laser.group.rotation, 'y', -m, m).listen();
	this.lrfolder.add(this.renderer.laser.group.rotation, 'z', -m, m).listen();
	this.gui.remember(this.renderer.laser.position);
	this.gui.remember(this.renderer.laser.rotation);
	this.lpfolder.close();
	this.lrfolder.close();
	*/
}

GUI.prototype.emissionSettings = function()
{
	this.emissionFolder = this.gui.addFolder('Emission');
	var renderer = this.renderer;
	var laser = renderer.laser;
	this.emissionSettings = {};
	settings = this.emissionSettings;
	settings.showLaserPointer = true;
	settings.spectrum = 'monochromatic';

	this.emissionFolder.add(settings, 'showLaserPointer').onChange( function(value) { laser.toggleVisibility(value); } );
	this.emissionFolder.add(laser, 'emissionRadius', 0.0, 10.0).onChange( function(value) { laser.setEmissionRadius(value);      renderer.reset(); } );
	this.emissionFolder.add(laser, 'emissionSpread', 0.0, 90.0).onChange( function(value) { laser.setEmissionSpreadAngle(value); renderer.reset(); } );
	this.emissionFolder.add(laser, 'emissionPower', 0.0, 10.0).onChange( function(value) { laser.setEmissionPower(value);       renderer.reset(); } );
	this.gui.remember(laser);

	// Spectrum selection
	var GUI = this;
	var spectrumObj = renderer.getLoadedSpectrum();
	var spectrumName = spectrumObj.getName();
	var spectra = renderer.getSpectra();
	var spectrumNames = Object.keys(spectra);

	settings["spectrum selection"] = spectrumName;
	this.emissionFolder.add(settings, 'spectrum selection', spectrumNames).onChange( function(spectrumName) {

						// remove gui for current spectrum
						var spectrumObj = renderer.getLoadedSpectrum();
						spectrumObj.eraseGui(GUI.emissionFolder);

						// load new scene
				 		renderer.loadSpectrum(spectrumName);

				 		// init gui for new scene
				 		spectrumObj = renderer.getLoadedSpectrum();
				 		spectrumObj.initGui(GUI.emissionFolder);
				 		
				 	} );

	spectrumObj.initGui(this.emissionFolder);

	this.gui.remember(settings);
	this.emissionFolder.open();
	return this.emissionFolder;
}


GUI.prototype.sceneSettings = function()
{
	this.sceneFolder = this.gui.addFolder('Scene');
	var renderer = this.renderer;
	var sceneObj = renderer.getLoadedScene();
	var sceneName = sceneObj.getName();
	scenes = renderer.getScenes();
	var sceneNames = Object.keys(scenes);

	// Scene selection menu
	settings = {};
	settings["scene selection"] = sceneName;
	var GUI = this;

	this.sceneFolder.add(settings, 'scene selection', sceneNames).onChange( function(sceneName) {

						// remove gui for current scene
						var sceneObj = renderer.getLoadedScene();
						sceneObj.eraseGui(GUI.sceneFolder);

						// load new scene
				 		renderer.loadScene(sceneName);

				 		// init gui for new scene
				 		sceneObj = renderer.getLoadedScene();
				 		sceneObj.initGui(GUI.sceneFolder);
				 		
				 	} );

	sceneObj.initGui(this.sceneFolder);
	this.sceneFolder.open();

	this.gui.remember(settings);
	return this.sceneFolder;
}


GUI.prototype.materialSettings = function()
{
	this.materialFolder = this.gui.addFolder('Material');
	var renderer = this.renderer;
	var materialObj = renderer.getLoadedMaterial();
	var materialName = materialObj.getName();
	materials = renderer.getMaterials();
	var materialNames = Object.keys(materials);

	// Material selection menu
	settings = {};
	settings["material selection"] = materialName;
	var GUI = this;

	this.materialFolder.add(settings, 'material selection', materialNames).onChange( function(materialName) {

						// remove gui for current material
						var materialObj = renderer.getLoadedMaterial();
						materialObj.eraseGui(GUI.materialFolder);

						// load new material
				 		renderer.loadMaterial(materialName);

				 		// init gui for new material
				 		materialObj = renderer.getLoadedMaterial();
				 		materialObj.initGui(GUI.materialFolder);
				 		
				 	} );
	
	materialObj.initGui(this.materialFolder);
	this.materialFolder.open();

	this.gui.remember(settings);
	return this.materialFolder;
}




