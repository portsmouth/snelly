

var GUI = function() 
{
	// Create dat gui
	this.gui = new dat.GUI();
	this.gui.domElement.id = 'gui';
	var gui = this.gui;
	this.visible = true;
	
	this.createSceneSettings();
	this.createMaterialSettings();
	this.createLightingSettings();
	this.createPathtracerSettings();
}

function updateDisplay(gui) 
{
    for (var i in gui.__controllers) {
        gui.__controllers[i].updateDisplay();
    }
    for (var f in gui.__folders) {
        updateDisplay(gui.__folders[f]);
    }
}

GUI.prototype.sync = function()
{
	updateDisplay(this.gui);
}

function hexToRgb(hex) 
{
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

	// @todo: add a basic AO and normals mode as well, useful for scene debugging.
	//var renderModes = ['normals', 'blinn'];
	//this.pathtracerFolder.add(pathtracer, 'renderMode', renderModes).onChange( function(renderMode) { pathtracer.reset(); });
	
	this.pathtracerFolder.add(pathtracer, 'exposure', 0.0, 50.0);
	this.pathtracerFolder.add(pathtracer, 'gamma', 0.0, 3.0);
	this.pathtracerFolder.add(pathtracer, 'whitepoint', 0.0, 2.0);
	this.pathtracerFolder.add(pathtracer, 'maxBounces', 1, 10).onChange( function(value) { pathtracer.maxBounces = Math.floor(value); pathtracer.reset(); } );
	this.pathtracerFolder.add(pathtracer, 'maxMarchSteps', 1, 1024).onChange( function(value) { pathtracer.maxMarchSteps = Math.floor(value); pathtracer.reset(); } );

	this.gui.remember(this.pathtracerSettings);
	this.pathtracerFolder.open();
}

GUI.prototype.createLightingSettings = function()
{
	this.emissionFolder = this.gui.addFolder('Lighting');
	var pathtracer = snelly.getPathtracer();
	this.emissionSettings = {};
	this.emissionSettings.spectrum = 'monochromatic';
	
	var skyPowerItem = this.emissionFolder.add(pathtracer, 'skyPower', 0.0, 10.0);
	skyPowerItem.onChange( function(value) 
	{ 
		snelly.camControls.enabled = false;
		var no_recompile = true;
		pathtracer.reset(no_recompile);
	} );
	skyPowerItem.onFinishChange( function(value) { snelly.camControls.enabled = true; } );

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
}

GUI.prototype.createSceneSettings = function()
{
	this.sceneFolder = this.gui.addFolder('Scene');
	var sceneObj = snelly.getScene();
	sceneObj.initGui(this);
	this.sceneFolder.open();
	this.gui.remember(this.sceneSettings);
}

GUI.prototype.addParameter = function(parameters, param)
{
	var name = param.name;
	var min  = param.min;
	var max  = param.max;
	var step = param.step;
	var recompile = param.recompile;
	var no_recompile = true;
	if (!(recompile==null || recompile==undefined)) no_recompile = !recompile;
	var item;
	if (step==null || step==undefined) { item = this.sceneFolder.add(parameters, name, min, max, step); }
	else                               { item = this.sceneFolder.add(parameters, name, min, max);       }
	item.onChange( function(value) { snelly.reset(no_recompile); snelly.camControls.enabled = false; } );
	item.onFinishChange( function(value) { snelly.camControls.enabled = true; } );
}

GUI.prototype.getSceneFolder = function()
{
	return this.sceneFolder;
}

GUI.prototype.createMaterialSettings = function()
{
	var GUI = this;

	var sceneObj = snelly.getScene();
	var sdfCode = sceneObj.sdf();

	// Metal settings
	if (sdfCode.indexOf("SDF_METAL(") !== -1)
	{
		this.metalFolder = this.gui.addFolder('Metal material');
		var metalObj = snelly.getLoadedMetal();
		var metalName = metalObj.getName();
		metals = snelly.getMetals();
		var metalNames = Object.keys(metals);
		this.metalMaterialSettings = {};
		this.metalMaterialSettings["metal"] = metalName;
		var metalItem = this.metalFolder.add(this.metalMaterialSettings, 'metal', metalNames);
		metalItem.onChange( function(materialName) {
							var materialObj = snelly.getLoadedMetal(); // remove gui for current material
							materialObj.eraseGui(GUI.metalFolder);
					 		snelly.loadMetal(materialName); // load new material
					 		materialObj = snelly.getLoadedMetal(); // init gui for new material
					 		materialObj.initGui(GUI.metalFolder);
				 	} );
		metalObj.initGui(this.metalFolder);
		this.metalFolder.open();
	}

	// Dielectric settings
	if (sdfCode.indexOf("SDF_DIELECTRIC(") !== -1)
	{
		this.dielectricFolder = this.gui.addFolder('Dielectric material');
		var dielectricObj = snelly.getLoadedDielectric();
		var dielectricName = dielectricObj.getName();
		dielectrics = snelly.getDielectrics();
		var dielectricNames = Object.keys(dielectrics);
		this.dielMaterialSettings = {};
		this.dielMaterialSettings["dielectric"] = dielectricName;
		var dielItem = this.dielectricFolder.add(this.dielMaterialSettings, 'dielectric', dielectricNames);
		dielItem.onChange( function(materialName) {
							var materialObj = snelly.getLoadedDielectric(); // remove gui for current material
							materialObj.eraseGui(GUI.dielectricFolder);
					 		snelly.loadDielectric(materialName); // load new material
					 		materialObj = snelly.getLoadedDielectric(); // init gui for new material
					 		materialObj.initGui(GUI.dielectricFolder);
					 	} );
		dielectricObj.initGui(this.dielectricFolder);
		this.dielectricFolder.open();
	}

	// Surface settings
	if (sdfCode.indexOf("SDF_SURFACE(") !== -1)
	{
		this.surfaceFolder = this.gui.addFolder('Surface material');
		var pathtracer = snelly.getPathtracer();

		this.surfaceFolder.diffuseReflectance = [pathtracer.surfaceDiffuseAlbedoRGB[0]*255.0, pathtracer.surfaceDiffuseAlbedoRGB[1]*255.0, pathtracer.surfaceDiffuseAlbedoRGB[2]*255.0];
		var diffItem = this.surfaceFolder.addColor(this.surfaceFolder, 'diffuseReflectance');
		diffItem.onChange( function(albedo) {
								if (typeof albedo==='string' || albedo instanceof String)
								{
									var color = hexToRgb(albedo);
									pathtracer.surfaceDiffuseAlbedoRGB[0] = color.r / 255.0;
									pathtracer.surfaceDiffuseAlbedoRGB[1] = color.g / 255.0;
									pathtracer.surfaceDiffuseAlbedoRGB[2] = color.b / 255.0;
								}
								else
								{
									pathtracer.surfaceDiffuseAlbedoRGB[0] = albedo[0] / 255.0;
									pathtracer.surfaceDiffuseAlbedoRGB[1] = albedo[1] / 255.0;
									pathtracer.surfaceDiffuseAlbedoRGB[2] = albedo[2] / 255.0;
								}
								pathtracer.surfaceDiffuseAlbedoXYZ = rgbToXyz(pathtracer.surfaceDiffuseAlbedoRGB);
								snelly.reset(true);
							} );

		this.surfaceFolder.specularReflectance = [pathtracer.surfaceSpecAlbedoRGB[0]*255.0, pathtracer.surfaceSpecAlbedoRGB[1]*255.0, pathtracer.surfaceSpecAlbedoRGB[2]*255.0];
		var specItem = this.surfaceFolder.addColor(this.surfaceFolder, 'specularReflectance');
		specItem.onChange( function(albedo) {
								if (typeof albedo==='string' || albedo instanceof String)
								{
									var color = hexToRgb(albedo);
									pathtracer.surfaceSpecAlbedoRGB[0] = color.r / 255.0;
									pathtracer.surfaceSpecAlbedoRGB[1] = color.g / 255.0;
									pathtracer.surfaceSpecAlbedoRGB[2] = color.b / 255.0;
								}
								else
								{
									pathtracer.surfaceSpecAlbedoRGB[0] = albedo[0] / 255.0;
									pathtracer.surfaceSpecAlbedoRGB[1] = albedo[1] / 255.0;
									pathtracer.surfaceSpecAlbedoRGB[2] = albedo[2] / 255.0;
								}
								pathtracer.surfaceSpecAlbedoXYZ = rgbToXyz(pathtracer.surfaceSpecAlbedoRGB);
								snelly.reset(true);
							} );

		this.roughness = pathtracer.surfaceRoughness;
		this.roughnessItem = this.surfaceFolder.add(this, 'roughness', 0.0, 0.1);
		this.roughnessItem.onChange( function(value) { pathtracer.surfaceRoughness = value; snelly.camControls.enabled = false; snelly.reset(true); } );
		this.roughnessItem.onFinishChange( function(value) { snelly.camControls.enabled = true; } );

		this.ior = pathtracer.surfaceIor;
		this.iorItem = this.surfaceFolder.add(this, 'ior', 0.0, 10.0);
		this.iorItem.onChange( function(value) { pathtracer.surfaceIor = value; snelly.camControls.enabled = false; snelly.reset(true); } );
		this.iorItem.onFinishChange( function(value) { snelly.camControls.enabled = true; } );

		this.surfaceFolder.open();
	}

	this.gui.remember(this.materialSettings);
}




