
function OceanScene(name, desc) 
{
	Scene.call(this, name, desc);

	// defaults
	this._settings.SEA_FREQ = 0.16;
	this._settings.SEA_HEIGHT = 0.6;
	this._settings.SEA_CHOPPY = 4.0;
	this._settings.SEA_TIME = 0.0;
}

// NB, every function is mandatory and must be defined.

OceanScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
OceanScene.prototype.sdf = function()
{
	return `
			// From https://www.shadertoy.com/view/Ms2SD1 by TDM
			uniform float SEA_FREQ;
			uniform float SEA_HEIGHT;
			uniform float SEA_CHOPPY;
			uniform float SEA_TIME;

			mat2 octave_m = mat2(1.6,1.2,-1.2,1.6);

			float hash( vec2 p ) 
			{
				float h = dot(p,vec2(127.1,311.7));	
			    return fract(sin(h)*43758.5453123);
			}

			float noise( in vec2 p ) 
			{
			    vec2 i = floor( p );
			    vec2 f = fract( p );	
				vec2 u = f*f*(3.0-2.0*f);
			    return -1.0+2.0*mix( mix( hash( i + vec2(0.0,0.0) ), 
			                     hash( i + vec2(1.0,0.0) ), u.x),
			                mix( hash( i + vec2(0.0,1.0) ), 
			                     hash( i + vec2(1.0,1.0) ), u.x), u.y);
			}

			float sea_octave(vec2 uv, float choppy) 
			{
			    uv += noise(uv);        
			    vec2 wv = 1.0-abs(sin(uv));
			    vec2 swv = abs(cos(uv));    
			    wv = mix(wv,swv,wv);
			    return pow(1.0-pow(wv.x * wv.y,0.65),choppy);
			}

			float ocean(vec3 p) 
			{
			    float freq = SEA_FREQ;
			    float amp = SEA_HEIGHT;
			    float choppy = SEA_CHOPPY;
			    vec2 uv = p.xz; uv.x *= 0.75;
			    float d, h = 0.0;    
			    const int ITER_GEOMETRY = 3;
			    for(int i = 0; i < ITER_GEOMETRY; i++) 
			    {        
			    	d = sea_octave((uv+SEA_TIME)*freq,choppy);
			    	d += sea_octave((uv-SEA_TIME)*freq,choppy);
			        h += d * amp;        
			    	uv *= octave_m; freq *= 1.9; amp *= 0.22;
			        choppy = mix(choppy,1.0,0.2);
			    }
			    return p.y - h;
			}

			float SDF(vec3 X)                     
			{                 
				return ocean(X);
			}                                     
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
OceanScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("SEA_FREQ", this._settings.SEA_FREQ);
	traceProgram.uniformF("SEA_HEIGHT", this._settings.SEA_HEIGHT);
	traceProgram.uniformF("SEA_CHOPPY", this._settings.SEA_CHOPPY);
	traceProgram.uniformF("SEA_TIME", this._settings.SEA_TIME);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
OceanScene.prototype.getScale = function()
{
	return 0.1;
}


// Initial cam position default for this scene
OceanScene.prototype.setCam = function(controls, camera)
{
	camera.position.set(10.0, 10.0, 10.0)
	controls.target.set(0.0, 0.0, 0.0);
}


// Initial laser position and direction defaults for this scene
OceanScene.prototype.setLaser = function(laser)
{
	laser.setPosition(new THREE.Vector3(0.0, 30.0, 0.0));
	laser.setDirection(new THREE.Vector3(0.0, -1.0, 0.0));
	Scene.prototype.setLaser.call(this, laser);
}


// set up gui and callbacks for this scene
OceanScene.prototype.initGui = function(parentFolder)
{
	this.itemSEA_FREQ = parentFolder.add(this._settings, 'SEA_FREQ', 0.0, 1.0);
	this.itemSEA_HEIGHT = parentFolder.add(this._settings, 'SEA_HEIGHT', 0.0, 2.0);
	this.itemSEA_CHOPPY = parentFolder.add(this._settings, 'SEA_CHOPPY', 0.0, 10.0);
	this.itemSEA_TIME = parentFolder.add(this._settings, 'SEA_TIME', 0.0, 10.0);

	this.itemSEA_FREQ.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_HEIGHT.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_CHOPPY.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_TIME.onChange( function(value) { snelly.reset(); } );
}

OceanScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.itemSEA_FREQ);
	parentFolder.remove(this.itemSEA_HEIGHT);
	parentFolder.remove(this.itemSEA_CHOPPY);
	parentFolder.remove(this.itemSEA_TIME);	
}










