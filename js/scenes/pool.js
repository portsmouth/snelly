
function PoolScene(name, desc) 
{
	Scene.call(this, name, desc);

	this._settings.SEA_FREQ = 0.16;
	this._settings.SEA_HEIGHT = 0.6;
	this._settings.SEA_DEPTH = 0.5;
	this._settings.SEA_CHOPPY = 4.0;
	this._settings.SEA_TIME = 0.0;
	this._settings.SEA_OCTAVES = 4;
}

// NB, every function is mandatory and must be defined.

PoolScene.prototype = Object.create(Scene.prototype);


// This defines a solid body, whose interior
// defined by the points with SDF<0.0, with a constant refractive index.
PoolScene.prototype.sdf = function()
{
	return `
			// Borrowed from https://www.shadertoy.com/view/Ms2SD1 by TDM
			uniform float SEA_FREQ;
			uniform float SEA_HEIGHT;
			uniform float SEA_DEPTH;
			uniform float SEA_CHOPPY;
			uniform float SEA_TIME;
			uniform float SEA_BOUNDS;

			float sdBox(vec3 X, vec3 bmin, vec3 bmax)                     
			{                            
				vec3 center = 0.5*(bmin + bmax);
				vec3 halfExtents = 0.5*(bmax - bmin);         
				vec3 d = abs(X-center) - halfExtents;
				return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));     
			} 

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
			    return pow(1.0-pow(wv.x * wv.y,0.65), choppy);
			}

			float ocean(vec3 p) 
			{
			    float freq = SEA_FREQ;
			    float amp = SEA_HEIGHT;
			    float choppy = SEA_CHOPPY;
			    vec2 uv = p.xz; uv.x *= 0.75;
			    float d, h = 0.0;    
			    for(int i = 0; i < ${Math.floor(this._settings.SEA_OCTAVES)}; i++) 
			    {        
			    	d  = sea_octave((uv+SEA_TIME)*freq,choppy);
			    	d += sea_octave((uv-SEA_TIME)*freq,choppy);
			        h += d * amp;        
			    	uv *= octave_m; freq *= 1.9; amp *= 0.22;
			        choppy = mix(choppy,1.0,0.2);
			    }
			    return p.y - h;
			}

			float SDF_DIELE(vec3 X)                     
			{                 
				vec3 bmax = vec3(SEA_BOUNDS,  4.0*SEA_HEIGHT, SEA_BOUNDS);
				vec3 bmin = vec3(-SEA_BOUNDS,     -SEA_DEPTH, -SEA_BOUNDS);	
				float bounds = sdBox(X, bmin, bmax);
				return opI(bounds, ocean(X));
			}                                  
	`;
}


// Called whenever this scene UI was switched to, or changed while active,
// and syncs the params of the trace shader to the current UI settings
PoolScene.prototype.syncShader = function(traceProgram)
{
	// (The shader parameter names here must be consistent with the GLSL sdf code defined above)
	traceProgram.uniformF("SEA_FREQ", this._settings.SEA_FREQ);
	traceProgram.uniformF("SEA_HEIGHT", this._settings.SEA_HEIGHT);
	traceProgram.uniformF("SEA_DEPTH", this._settings.SEA_DEPTH);
	traceProgram.uniformF("SEA_CHOPPY", this._settings.SEA_CHOPPY);
	traceProgram.uniformF("SEA_TIME", this._settings.SEA_TIME);

	var b = 50.0;
	traceProgram.uniformF("SEA_BOUNDS", b);
}

// Gives the raytracer some indication of (rough) scene size, so it
// can set tolerances appropriately.
PoolScene.prototype.getScale = function()
{
	return 50.0;
}


PoolScene.prototype.getBox = function()
{
	var h = this._settings.SEA_HEIGHT;
	var d = this._settings.SEA_DEPTH;
	var n = this._settings.SEA_OCTAVES;
	var b = 50.0;
	var min = new THREE.Vector3(-b,   -d, -b);
	var max = new THREE.Vector3( b,  h*n,  b);
	return new THREE.Box3(min, max);
}



// Initial cam position default for this scene
PoolScene.prototype.init = function(controls, camera, laser)
{
	laser.setPosition(new THREE.Vector3(0.00000, 10.0000, 0.00000));
	laser.setDirection(new THREE.Vector3(1.22465e-16, -1.00000, 0.00000));
	laser.setEmissionRadius(3.00000);
	laser.setEmissionSpreadAngle(45.0000);
	controls.target.set(5.55635, -1.05045, 1.96451);
	camera.position.set(50.8994, 20.0778, 35.8554);
}


// set up gui and callbacks for this scene
PoolScene.prototype.initGui = function(parentFolder)
{
	this.itemSEA_FREQ = parentFolder.add(this._settings, 'SEA_FREQ', 0.0, 0.5);
	this.itemSEA_HEIGHT = parentFolder.add(this._settings, 'SEA_HEIGHT', 0.0, 3.0);
	this.itemSEA_DEPTH = parentFolder.add(this._settings, 'SEA_DEPTH', 0.0, 50.0);
	this.itemSEA_CHOPPY = parentFolder.add(this._settings, 'SEA_CHOPPY', 0.0, 10.0);
	this.itemSEA_TIME = parentFolder.add(this._settings, 'SEA_TIME', 0.0, 30.0);
	this.itemSEA_OCTAVES = parentFolder.add(this._settings, 'SEA_OCTAVES', 1, 5, 1);

	this.itemSEA_FREQ.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_HEIGHT.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_DEPTH.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_CHOPPY.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_TIME.onChange( function(value) { snelly.reset(); } );
	this.itemSEA_OCTAVES.onChange( function(value) { snelly.reset(); } );
}

PoolScene.prototype.eraseGui = function(parentFolder)
{
	parentFolder.remove(this.itemSEA_FREQ);
	parentFolder.remove(this.itemSEA_HEIGHT);
	parentFolder.remove(this.itemSEA_DEPTH);
	parentFolder.remove(this.itemSEA_CHOPPY);
	parentFolder.remove(this.itemSEA_TIME);	
	parentFolder.remove(this.itemSEA_OCTAVES);	
}










