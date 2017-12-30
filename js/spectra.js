

////////////////////////////////////////////////////////
// Spectrum
////////////////////////////////////////////////////////

function Spectrum(name, desc)
{
    this._name = name;
    this._desc = desc;
}

Spectrum.prototype.getName = function()
{
    return this._name;
}

Spectrum.prototype.getDesc = function()
{
    return this._desc;
}

Spectrum.prototype.spectrum = function(wavelength)
{
    console.log('Spectrum.prototype.spectrum, error');
}

Spectrum.prototype.inverseCDF = function(minwavelength, maxwavelength, numsamples)
{
    // Compute CDF of spectrum
    var cdf  = new Float32Array(numsamples+1);
    var dl = (maxwavelength - minwavelength)/numsamples;
    cdf[0] = 0.0;
    for (var n=0; n<numsamples; ++n)
    {
        var wavelength = minwavelength + dl*n;
        cdf[n+1] = cdf[n] + this.spectrum(wavelength); 
    }

    // Normalize CDF
    var cdfSum = cdf[numsamples];
    for (var n=0; n<numsamples; ++n) cdf[n+1] /= cdfSum;
    cdf[numsamples] = 1.0;

    // Numerically invert CDF
    var icdf = new Float32Array(4*numsamples);
    var cdfIdx = 0;
    for (var n=0; n<4*numsamples; ++n)
    {
        var xi = Math.min((n+1)/(4*numsamples), 1.0);
        while (cdf[cdfIdx] < xi) cdfIdx++;
        cdfIdx = Math.max(1, cdfIdx);
        var xiLo = cdf[cdfIdx-1];
        var xiHi = cdf[cdfIdx];
        var icdfLo = (cdfIdx-1.0)/numsamples;
        var icdfHi =     (cdfIdx)/numsamples;
        var icdfVal = icdfLo + (1.0/numsamples)*(xi - xiLo)/(xiHi - xiLo);
        icdf[n] = icdfVal;
    }
    return icdf;
}

////////////////////////////////////////////////////////
// Monochromatic
////////////////////////////////////////////////////////

function MonochromaticSpectrum(name, desc, wavelength)
{
    Spectrum.call(this, name, desc);
    this.wavelength = wavelength;
}

MonochromaticSpectrum.prototype = Object.create(Spectrum.prototype);

// A delta function spectrum has specialized inverse CDF generation:
MonochromaticSpectrum.prototype.inverseCDF = function(minwavelength, maxwavelength, numsamples)
{
    var wavelength = Math.max(minwavelength, Math.min(this.wavelength, maxwavelength));
    var spectrumOffset = (wavelength - minwavelength) / (maxwavelength - minwavelength);
    var icdf = new Float32Array(4*numsamples);
    for (var n=0; n<4*numsamples; ++n)
    {
        icdf[n] = spectrumOffset;
    }
    return icdf;
}

// set up gui and callbacks for this spectrum
MonochromaticSpectrum.prototype.initGui = function(parentFolder)
{
    ME = this;
    this.wavelengthItem = parentFolder.add(this, 'wavelength', 390.0, 750.0);
    this.wavelengthItem.onChange( function(value) 
    { 
        snelly.controls.enabled = false;
        var no_recompile = true;
        snelly.loadSpectrum(ME.getName());
        snelly.reset(no_recompile);
    } );
    this.wavelengthItem.onFinishChange( function(value) { snelly.controls.enabled = true; } );
}

MonochromaticSpectrum.prototype.eraseGui = function(parentFolder)
{
    parentFolder.remove(this.wavelengthItem);
}



////////////////////////////////////////////////////////
// Flat
////////////////////////////////////////////////////////

function FlatSpectrum(name, desc, bandmin, bandmax)
{
    Spectrum.call(this, name, desc);
    this['Minimum wavelength'] = bandmin;
    this['Maximum wavelength'] = bandmax;
}

FlatSpectrum.prototype = Object.create(Spectrum.prototype);

FlatSpectrum.prototype.spectrum = function(wavelength)
{
    if (wavelength<this['Minimum wavelength']) return 0.0;
    if (wavelength>this['Maximum wavelength']) return 0.0;
    return 1.0;
}

// set up gui and callbacks for this spectrum
FlatSpectrum.prototype.initGui = function(parentFolder)
{
    ME = this;
    var dw = (750.0-390.0)/1024.0;
    this.minItem = parentFolder.add(this, 'Minimum wavelength', 390.0, 750.0, 0.1);
    this.minItem.onChange( function(value) 
    { 
        snelly.camControls.enabled = false;
        if (ME['Minimum wavelength'] > ME['Maximum wavelength']-dw) 
            ME['Minimum wavelength'] = ME['Maximum wavelength']-dw;
        var no_recompile = true;
        snelly.loadSpectrum(ME.getName());
        snelly.reset(no_recompile);
    });
    this.minItem.onFinishChange( function(value) { snelly.camControls.enabled = true; } );

    this.maxItem = parentFolder.add(this, 'Maximum wavelength', 390.0, 750.0, 0.1);
    this.maxItem.onChange( function(value) 
    {
        snelly.camControls.enabled = false; 
        if (ME['Maximum wavelength'] < ME['Minimum wavelength']+dw) 
            ME['Maximum wavelength'] = ME['Minimum wavelength']+dw;
        var no_recompile = true;
        snelly.loadSpectrum(ME.getName());
        snelly.reset(no_recompile);
    });
    this.maxItem.onFinishChange( function(value) { snelly.camControls.enabled = true; } );
}

FlatSpectrum.prototype.eraseGui = function(parentFolder)
{
    parentFolder.remove(this.minItem);
    parentFolder.remove(this.maxItem);
}


////////////////////////////////////////////////////////
// Blackbody
////////////////////////////////////////////////////////

function BlackbodySpectrum(name, desc, temperature) // temperature in Kelvin
{
    Spectrum.call(this, name, desc);
    this.skyTemperature = temperature;
}

BlackbodySpectrum.prototype = Object.create(Spectrum.prototype);

BlackbodySpectrum.prototype.spectrum = function(wavelength)  // wavelength in nm
{
    var boltzmann_factor = 1.43877737467e7 / (wavelength*this.skyTemperature);
    var l = wavelength/360.0; // wavelength relative to 360nm
    return 1.0 / (l*l*l*l*l*(Math.exp(boltzmann_factor) - 1.0));
}

// set up gui and callbacks for this spectrum
BlackbodySpectrum.prototype.initGui = function(parentFolder)
{
    ME = this;
    this.temperatureItem = parentFolder.add(this, 'skyTemperature', 300.0, 15000.0);
    this.temperatureItem.onChange( function(value) 
    { 
        snelly.camControls.enabled = false;
        this.skyTemperature = value;
        snelly.loadSpectrum(ME.getName());
        snelly.reset(true);
    } );
    this.temperatureItem.onFinishChange( function(value) { snelly.camControls.enabled = true; } );
}

BlackbodySpectrum.prototype.eraseGui = function(parentFolder)
{
    parentFolder.remove(this.temperatureItem);
}



