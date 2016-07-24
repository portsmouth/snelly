

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
	folder = this.gui.addFolder('Renderer');

	var renderer = this.renderer;
	renderer.renderSettings = {};
	settings = renderer.renderSettings;

	settings.exposure = 100.0;
	settings.maxPathLength = 32;
	settings.maxMarchSteps = 512;
	settings.photonsPerFrame = renderer.raySize;
	settings.showSurface = false;
	settings.surfaceAlpha = 0.5;
	
	folder.add(settings, 'exposure', 0.0, 1000.0, 0.01);
	folder.add(settings, 'maxPathLength', 1, 1024).onChange( function(value) { renderer.renderSettings.maxPathLength = Math.floor(value); renderer.reset(); } );
	folder.add(settings, 'maxMarchSteps', 1, 1024).onChange( function(value) { renderer.renderSettings.maxMarchSteps = Math.floor(value); renderer.reset(); } );
	folder.add(settings, 'photonsPerFrame', 1, 1024).onChange( function(value) { renderer.raySize = Math.floor(value); renderer.initRayStates(); } );
	folder.add(settings, 'showSurface');
	folder.add(settings, 'surfaceAlpha');
	this.gui.remember(settings);

	folder.open();
	return folder;
}

GUI.prototype.registerCamera = function()
{
	rs = this.rs;
	var m = 100.0;
	this.cpfolder = rs.addFolder('Camera Position');
	this.cpfolder.add(this.renderer.camera.position, 'x', -m, m).listen();
	this.cpfolder.add(this.renderer.camera.position, 'y', -m, m).listen();
	this.cpfolder.add(this.renderer.camera.position, 'z', -m, m).listen();
	this.crfolder = rs.addFolder('Camera Rotation');
	this.crfolder.add(this.renderer.camera.rotation, 'x', -m, m).listen();
	this.crfolder.add(this.renderer.camera.rotation, 'y', -m, m).listen();
	this.crfolder.add(this.renderer.camera.rotation, 'z', -m, m).listen();
	this.gui.remember(this.renderer.camera.position);
	this.gui.remember(this.renderer.camera.rotation);
	this.cpfolder.close();
	this.crfolder.close();
}

GUI.prototype.registerLaser = function()
{
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
}

GUI.prototype.emissionSettings = function()
{
	folder = this.gui.addFolder('Emission');

	var renderer = this.renderer;
	var laser = renderer.laser;

	this.emissionSettings = {};
	settings = this.emissionSettings;

	settings.radius = laser.getEmissionRadius();
	settings.spread = laser.getEmissionSpreadAngle();
	settings.showLaserPointer = true;
	settings.spectrum = 'monochromatic';

	folder.add(settings, 'showLaserPointer').onChange( function(value) { laser.toggleVisibility(value); } );
	folder.add(laser, 'emissionRadius', 0.0, 10.0).onChange( function(value) { laser.setEmissionRadius(value);      renderer.reset(); } );
	folder.add(laser, 'emissionSpread', 0.0, 90.0).onChange( function(value) { laser.setEmissionSpreadAngle(value); renderer.reset(); } );
	folder.add(settings, 'spectrum', ['monochromatic', 'color', 'flat', 'blackbody']);

	this.gui.remember(settings);

	spectrumParams = folder.addFolder('spectrumParams');
	spectrumParams.close();

		monochromaticFolder = spectrumParams.addFolder('monochromatic');
		settings.monochromaticSettings = {wavelength: 500.0};
		monochromaticFolder.add(settings.monochromaticSettings, 'wavelength', 450.0, 650.0);
		this.gui.remember(settings.monochromaticSettings);

		colorFolder = spectrumParams.addFolder('color');
		settings.colorSettings = {colorPicker: { h: 350, s: 0.9, v: 0.3 }};
		colorFolder.addColor(settings.colorSettings, 'colorPicker');
		this.gui.remember(settings.colorSettings);

		flatFolder = spectrumParams.addFolder('flat');
		settings.flatSettings = {minwavelength: 450.0, maxwavelength: 650.0};
		flatFolder.add(settings.flatSettings, 'minwavelength', 450.0, 650.0);
		flatFolder.add(settings.flatSettings, 'maxwavelength', 450.0, 650.0);
		this.gui.remember(settings.flatSettings);

		blackbodyFolder = spectrumParams.addFolder('blackbody');
		settings.blackbodySettings = {temperature: 6000.0};
		blackbodyFolder.add(settings.blackbodySettings, 'temperature', 1000.0, 15000.0);
		this.gui.remember(settings.blackbodySettings);

	folder.open();
	return folder;
}

GUI.prototype.sceneSettings = function()
{
	folder = this.gui.addFolder('Scene');

	var renderer = this.renderer;
	this.sceneSettings = {};
	settings = this.sceneSettings;

	var sceneObj = renderer.getLoadedScene();
	settings.scene = sceneObj.getName();

	scenes = renderer.getScenes();
	var sceneNames = Object.keys(scenes);

	folder.add(settings, 'scene', sceneNames).onChange( function(sceneName) { 
				 		console.log('Scene selection changed to ' + sceneName);
				 		renderer.loadScene(sceneName);
				 	} );

	this.gui.remember(settings);

	return folder;
}

GUI.prototype.materialSettings = function()
{
	folder = this.gui.addFolder('Material');

	var renderer = this.renderer;
	this.materialSettings = {};
	settings = this.materialSettings;

	var materialObj = renderer.getLoadedScene();
	settings.material = sceneObj.materialObj();
	
	materials = renderer.getMaterials();
	var materialNames = Object.keys(materials);

	folder.add(settings, 'material', materialNames).onChange( function(materialName) { 
				 		console.log('Material selection changed to ' + materialName);
				 		renderer.loadMaterial(materialName);
				 	} );

	this.gui.remember(settings);

	return folder;
}























