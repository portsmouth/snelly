
function WavesScene(name, desc) 
{
	Scene.call(this, name, desc);

	this._settings.waveFreq1 = 1.0;
	this._settings.waveFreq2 = 0.3;
	this._settings.waveHeight1 = 0.1;
	this._settings.waveHeight2 = 0.1;
	this._settings.length = 10.0;
	this._settings.width = 10.0;
	this._settings.angle1 = 0.0;
	this._settings.angle2 = 45.0;
	this._settings.phase1 = 0.0;
	this._settings.phase2 = 270.0;
	this._settings.depth = 3.0;
}

// NB, every function is mandatory and must be defined.

WavesScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
WavesScene.prototype.sdf = function()
{
	return `
			uniform float _waveFreq1;
			uniform float _waveFreq2;
			uniform float _waveHeight1;
			uniform float _waveHeight2;
			uniform float _angle1;
			uniform float _angle2;
			uniform float _phase1;
			uniform float _phase2;

			uniform float _length;
			uniform float _width;
			uniform float _depth;

			float ocean(vec3 p) 
			{
				float wx1 = 2.0*M_PI*_waveFreq1*cos(_angle1);
				float wz1 = 2.0*M_PI*_waveFreq1*sin(_angle1);
			    float h1 = _waveHeight1 * sin(wx1*p.x+wz1*p.z + _phase1);

			    float wx2 = 2.0*M_PI*_waveFreq2*cos(_angle2);
				float wz2 = 2.0*M_PI*_waveFreq2*sin(_angle2);
				float h2 = _waveHeight2 * sin(wx2*p.x+wz2*p.z + _phase2);

			    return p.y - (h1 + h2);
			}

			float SDF_DIELE(vec3 X)                     
			{                 
				vec3 bmax = vec3(_length,  4.0*(_waveHeight1 + _waveHeight2), _width);
				vec3 bmin = vec3(-_length,     -_depth, -_width);	
				float bounds = sdBox(X, bmin, bmax);
				return opI(bounds, ocean(X));
			}    
				
			float SDF_METAL(vec3 X) { return HUGE_VAL; }
			float SDF_DIFFU(vec3 X) { return sdBox(X, vec3(-100.0, -2.5, -100.0), vec3(100.0, -2.0, 100.0)); }                                
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
WavesScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("_waveFreq1", this._settings.waveFreq1);
	traceProgram.uniformF("_waveFreq2", this._settings.waveFreq2);
	traceProgram.uniformF("_waveHeight1", this._settings.waveHeight1);
	traceProgram.uniformF("_waveHeight2", this._settings.waveHeight2);
	traceProgram.uniformF("_angle1", this._settings.angle1*Math.PI/180.0);
	traceProgram.uniformF("_angle2", this._settings.angle2*Math.PI/180.0);
	traceProgram.uniformF("_phase1", this._settings.phase1*Math.PI/180.0);
	traceProgram.uniformF("_phase2", this._settings.phase2*Math.PI/180.0);

	traceProgram.uniformF("_length", this._settings.length);
	traceProgram.uniformF("_width", this._settings.width);
	traceProgram.uniformF("_depth", this._settings.depth);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
WavesScene.prototype.getScale = function()
{
	return Math.max(this._settings.length, this._settings.width, this._settings.depth);
}

/*
WavesScene.prototype.getBox = function()
{

}
*/


// Initial cam position default for this scene
WavesScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(-7.01068, 2.49419, 1.87921));
	laser.setTarget(new THREE.Vector3(-2.65843, 0.0817634, 6.21013));
	laser.setEmissionRadius(0.0100000);
	laser.setEmissionSpreadAngle(0.00000);
	controls.target.set(2.11526, -4.62020, -2.34717);
	camera.position.set(1.21316, 4.38071, 20.1111);
}


// set up gui and callbacks for this scene
WavesScene.prototype.initGui = function(parentFolder)
{
	this.itemFreq1 = parentFolder.add(this._settings, 'waveFreq1', 0.0, 2.0);
	this.itemHeight1 = parentFolder.add(this._settings, 'waveHeight1', 0.0, 0.3);
	this.itemAngle1 = parentFolder.add(this._settings, 'angle1', 0.0, 90.0);
	this.itemPhase1 = parentFolder.add(this._settings, 'phase1', 0.0, 360.0);

	this.itemFreq2 = parentFolder.add(this._settings, 'waveFreq2', 0.0, 2.0);
	this.itemHeight2 = parentFolder.add(this._settings, 'waveHeight2', 0.0, 0.3);
	this.itemAngle2 = parentFolder.add(this._settings, 'angle2', 0.0, 90.0);
	this.itemPhase2 = parentFolder.add(this._settings, 'phase2', 0.0, 360.0);

	this.itemLength = parentFolder.add(this._settings, 'length', 0.0, 30.0);
	this.itemWidth = parentFolder.add(this._settings, 'width', 0.0, 30.0);
	this.itemDepth = parentFolder.add(this._settings, 'depth', 0.0, 50.0);

	this.itemFreq1.onChange( function(value) { snelly.reset(); } );
	this.itemFreq2.onChange( function(value) { snelly.reset(); } );
	this.itemHeight1.onChange( function(value) { snelly.reset(); } );
	this.itemHeight2.onChange( function(value) { snelly.reset(); } );
	this.itemAngle1.onChange( function(value) { snelly.reset(); } );
	this.itemAngle2.onChange( function(value) { snelly.reset(); } );
	this.itemPhase1.onChange( function(value) { snelly.reset(); } );
	this.itemPhase2.onChange( function(value) { snelly.reset(); } );
	this.itemLength.onChange( function(value) { snelly.reset(); } );
	this.itemWidth.onChange( function(value) { snelly.reset(); } );
	this.itemDepth.onChange( function(value) { snelly.reset(); } );
}

WavesScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.itemFreq1);
	parentFolder.remove(this.itemFreq2);
	parentFolder.remove(this.itemHeight1);
	parentFolder.remove(this.itemHeight2);
	parentFolder.remove(this.itemAngle1);
	parentFolder.remove(this.itemAngle2);
	parentFolder.remove(this.itemPhase1);
	parentFolder.remove(this.itemPhase2);

	parentFolder.remove(this.itemLength);
	parentFolder.remove(this.itemWidth);
	parentFolder.remove(this.itemDepth);
}










