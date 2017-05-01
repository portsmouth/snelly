

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
	var pathtracer = snelly.getPathtracer();
	this.emissionSettings = {};
	this.emissionSettings.spectrum = 'monochromatic';
	
	this.emissionFolder.add(pathtracer, 'skyPower', 0.0, 1.0).onChange( function(value) 
	{ 
		laser.setSkyPower(value);
		pathtracer.reset();
	} );

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

	this.gui.remember(this.emissionSettings);
	this.emissionFolder.open();
	return this.emissionFolder;
}


GUI.prototype.createSceneSettings = function()
{
	this.sceneFolder = this.gui.addFolder('Scene');

	var sceneObj = snelly.getScene();
	var sceneName = sceneObj.getName();

	sceneObj.initGui(this.sceneFolder);
	this.sceneFolder.open();

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




