
var PathtracerState = function(width, height)
{
    var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count
    var rngData      = new Float32Array(width*height*4); // Random number seed
    for (var i=0; i<width*height; ++i)
    {
        rngData[i*4 + 0] = Math.random()*4194167.0;
        rngData[i*4 + 1] = Math.random()*4194167.0;
        rngData[i*4 + 2] = Math.random()*4194167.0;
        rngData[i*4 + 3] = Math.random()*4194167.0;
    }
    this.radianceTex = new GLU.Texture(width, height, 4, true, false, true, radianceData);
    this.rngTex      = new GLU.Texture(width, height, 4, true, false, true, rngData);
}

PathtracerState.prototype.bind = function(shader)
{
    this.radianceTex.bind(0);
    this.rngTex.bind(1);
    shader.uniformTexture("Radiance", this.radianceTex);
    shader.uniformTexture("RngData", this.rngTex);
}

PathtracerState.prototype.attach = function(fbo)
{
    var gl = GLU.gl;
    fbo.attachTexture(this.radianceTex, 0);
    fbo.attachTexture(this.rngTex, 1);
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE)
    {
        GLU.fail("Invalid framebuffer");
    }
}

PathtracerState.prototype.detach = function(fbo)
{
    var gl = GLU.gl;
    fbo.detachTexture(0);
    fbo.detachTexture(1);
}

PathtracerState.prototype.clear = function(fbo)
{
    // clear radiance buffer
    var gl = GLU.gl;
    fbo.bind();
    fbo.drawBuffers(1);
    fbo.attachTexture(this.radianceTex, 0);
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    fbo.unbind();
}

/**
* Interface to the renderer. The rendering modes available are:
*  - 'pt': pathtracer (uni-directional)
*  - 'ao': ambient occlusion, colored via {@link Surface} material diffuse albedo modulated by the `SURFACE_DIFFUSE_REFLECTANCE` shader function
*  - 'normals': view normal at first hit as a color
*  - 'firsthit': first-hit only, and {@link Surface} material only
* @constructor
* @property {number} width                     - (if not specified, fits to window)
* @property {number} height                    - (if not specified, fits to window)
* @property {String} [renderMode='pt']         - rendering mode (either 'pt', 'ao', 'normals', 'firsthit')
* @property {number} [maxSamplesPerFrame=1]    - maximum number of per-pixel samples per frame
* @property {number} [maxBounces=3]            - maximum number of surface bounces
* @property {number} [maxScatters=2]           - maximum number of volumetric scatterings
* @property {number} [maxMarchSteps=256]       - maximum number of raymarching steps per path segment
* @property {number} [maxVolumeSteps=256]      - maximum number of Woodcock tracking steps per path segment
* @property {number} [maxStepsIsMiss=true]     - whether rays which exceed max step count are considered hits or misses
* @property {number} [interactive=true]        - if enabled, tries to maintain interactive frame rate at the expense of more noise
* @property {number} [goalFPS=10.0]            - sampling will adjust to try to match goal FPS
* @property {number} [minsSPPToRedraw=0.0]     - if >0.0, renderer will not redraw until the specified SPP have been accumulated
* @property {number} [radianceClamp=3.0]       - clamp radiance to (10^) this max value, for firefly reduction
* @property {number} [wavelengthSamples=256]   - number of samples to take over visible wavelength range
* @property {number} [exposure=0.0]            - exposure, on a log scale
* @property {number} [gamma=2.2]               - display gamma correction
* @property {number} [contrast=1.0]            - tonemapping contrast
* @property {number} [saturation=1.0]          - tonemapping saturation
* @property {number} [skyPower=4.0]            - sky power (arbitrary units)
* @property {number} [skyTemperature=6000]     - sky temperature (in Kelvin)
* @property {number} [envMapRotation=0.0]      - env map rotation about pole in degrees (0 to 360)
* @property {number} [envMapVisible=true]      - whether env map is visible to primary rays (otherwise black)
* @property {number} [sunPower=1.0]            - sun power (arbitrary units)
* @property {Array}  [sunColor]                - sun color
* @property {number} [sunAngularSize=5.0]      - sun angular size (degrees)
* @property {number} [sunLatitude=50.0]        - sun latitude (degrees)
* @property {number} [sunLongitude=0.0]        - sun longitude (degrees)
* @property {number} [sunVisibleDirectly=true] - whether sun is directly visible
* @property {number} [shadowStrength=1.0]     - if <1.0, areas in shadow are not completely dark (provided mostly to allow rendering of occluded areas, e.g. fractals)
*/
var Renderer = function()
{
    this.gl = GLU.gl;
    var gl = GLU.gl;

    var render_canvas = snelly.render_canvas;
    render_canvas.width  = window.innerWidth;
    render_canvas.height = window.innerHeight;
    this._width = render_canvas.width;
    this._height = render_canvas.height;

    // Initialize pathtracing buffers and programs
    this.currentState = 0;
    this.pathStates = [new PathtracerState(this._width, this._height),
                       new PathtracerState(this._width, this._height)];
    this.fbo = null;
    this.aoProgram           = null;
    this.normalsProgram      = null;
    this.firsthitProgram     = null;
    this.pathtraceAllProgram = null;
    this.tonemapProgram      = null;
    this.pickProgram         = null;

    // Internal properties (@todo: use underscore to make this more explicit?)
    this.numSamples = 0;
    this.numFramesSinceReset = 0;
    this.numFramesSinceInit = 0;
    this.skipProbability = 0.0;
    this.frametime_measure_ms = 0.0;
    this.spp = 0.0;

    // Default user-adjustable properties
    this.renderMode = 'pt';
    this.maxSamplesPerFrame = 1;
    this.maxBounces = 3;
    this.maxScatters = 2;
    this.maxMarchSteps = 256;
    this.maxVolumeSteps = 256;
    this.radianceClamp = 3.0;
    this.wavelengthSamples = 256;
    this.skyPower = 1.0;
    this.skyTemperature = 6000.0;
    this.sunPower = 1.0;
    this.sunAngularSize = 5.0;
    this.sunColor = [1.0,0.8,0.5];
    this.sunLatitude = 50.0;
    this.sunLongitude = 0.0;
    this.sunVisibleDirectly = true;
    this.exposure = 0.0;
    this.gamma = 2.2;
    this.contrast = 1.0;
    this.saturation = 1.0;
    this.goalFPS = 20.0;
    this.minsSPPToRedraw = 0.0;
    this.envMapVisible = true;
    this.envMapRotation = 0.0;
    this.shadowStrength = 1.0;
    this.maxStepsIsMiss = true;
    this.interactive = true;

    // Load shaders
    this.shaderSources = GLU.resolveShaderSource({
        'pathtracer': {'v': 'pathtracer-vertex-shader', 'f': 'pathtracer-fragment-shader'},
        'ao':         {'v': 'ao-vertex-shader',         'f': 'ao-fragment-shader'        },
        'normals':    {'v': 'normals-vertex-shader',    'f': 'normals-fragment-shader'   },
        'firsthit':   {'v': 'firsthit-vertex-shader',   'f': 'firsthit-fragment-shader'  },
        'tonemapper': {'v': 'tonemapper-vertex-shader', 'f': 'tonemapper-fragment-shader'},
        'pick':       {'v': 'pick-vertex-shader',       'f': 'pick-fragment-shader'}
    });

    this.filterPrograms = null;

    // load env map
    this.loaded = true;
    var sceneObj = snelly.getScene();
    this.envMap = null;
    if (typeof(sceneObj.envMap) !== "undefined")
      {
          var url = sceneObj.envMap();
          if (typeof(url) != "undefined" && url != "")
          {
              var pathtracer = this;
              this.loaded = false;
              (function() { GLU.loadImageAndCreateTextureInfo(url,
                    function(imgInfo)
                    {
                        pathtracer.loaded =  true;
                        pathtracer.envMap = imgInfo;
                    });
              })(pathtracer.loaded);
          }
      }

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.blendFunc(gl.ONE, gl.ONE);

    this.quadVbo = this.createQuadVbo();
    this.fbo = new GLU.RenderTarget();

    // Trigger initial buffer generation
    this.resize(this._width, this._height);
}

Renderer.prototype.createQuadVbo = function()
{
    var vbo = new GLU.VertexBuffer();
    vbo.addAttribute("Position", 3, this.gl.FLOAT, false);
    vbo.addAttribute("TexCoord", 2, this.gl.FLOAT, false);
    vbo.init(4);
    vbo.copy(new Float32Array([
         1.0,  1.0, 0.0, 1.0, 1.0,
        -1.0,  1.0, 0.0, 0.0, 1.0,
        -1.0, -1.0, 0.0, 0.0, 0.0,
         1.0, -1.0, 0.0, 1.0, 0.0
    ]));
    return vbo;
}


/**
* Restart accumulating samples.
* @param {Boolean} [no_recompile=false] - set to true if shaders need recompilation too
*/
Renderer.prototype.reset = function(no_recompile = false)
{
    this.numSamples = 0;
    this.numFramesSinceReset = 0;
    if (!no_recompile) this.compileShaders();
    this.currentState = 0;
    this.pathStates[this.currentState].clear(this.fbo);
    this.pathStates[this.currentState+1].clear(this.fbo);
}

/**
* Read access to the current per-pixel sample count average.
* @returns {number} - the current sample count average
*/
Renderer.prototype.getSPP = function()
{
    return this.spp;
}

Renderer.prototype.compileShaders = function()
{
    // Inject code for the current scene SDF:
    var sceneObj = snelly.getScene();
    if (sceneObj == null) return;
    if (typeof sceneObj.shader == "undefined") { GLU.fail('Scene must define a "shader" function!'); }
    var shader = sceneObj.shader();

    // Insert default functions for those not implemented
    if (shader.indexOf("INIT(")                            == -1) { shader += `\n void INIT() {}\n`; }

    let hasSurface = true;
    let hasSurfaceNM = false;
    if (shader.indexOf("SDF_SURFACE(")                     == -1) { hasSurface = false; }
    if (shader.indexOf("SURFACE_DIFFUSE_REFLECTANCE(")     == -1) { shader += `\n vec3 SURFACE_DIFFUSE_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
    if (shader.indexOf("SURFACE_SPECULAR_REFLECTANCE(")    == -1) { shader += `\n vec3 SURFACE_SPECULAR_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
    if (shader.indexOf("SURFACE_ROUGHNESS(")               == -1) { shader += `\n float SURFACE_ROUGHNESS(in float roughness, in vec3 X, in vec3 N) { return roughness; }\n`; }
    if (shader.indexOf("SURFACE_NORMAL_MAP(")               > -1) { hasSurfaceNM = true; }

    let hasMetal = true;
    let hasMetalNM = false;
    if (shader.indexOf("SDF_METAL(")                       == -1) { hasMetal = false; }
    if (shader.indexOf("METAL_SPECULAR_REFLECTANCE(")      == -1) { shader += `\n vec3 METAL_SPECULAR_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
    if (shader.indexOf("METAL_ROUGHNESS(")                 == -1) { shader += `\n float METAL_ROUGHNESS(in float roughness, in vec3 X, in vec3 N) { return roughness; }\n`; }
    if (shader.indexOf("METAL_NORMAL_MAP(")                 > -1) { hasMetalNM = true; }

    let hasDielectric = true;
    let hasDielectricNM = false;
    if (shader.indexOf("SDF_DIELECTRIC(")                  == -1) { hasDielectric = false; }
    if (shader.indexOf("DIELECTRIC_SPECULAR_REFLECTANCE(") == -1) { shader += `\n vec3 DIELECTRIC_SPECULAR_REFLECTANCE(in vec3 C, in vec3 X, in vec3 N, in vec3 V) { return C; }\n`; }
    if (shader.indexOf("DIELECTRIC_ROUGHNESS(")            == -1) { shader += `\n float DIELECTRIC_ROUGHNESS(in float roughness, in vec3 X, in vec3 N) { return roughness; }\n`; }
    if (shader.indexOf("DIELECTRIC_NORMAL_MAP(")            > -1) { hasDielectricNM = true; }

    let hasVolume = true;
    if (shader.indexOf("SDF_VOLUME(")                      == -1) { hasVolume= false; }
    if (shader.indexOf("VOLUME_EXTINCTION(")               == -1) { shader += `\n float VOLUME_EXTINCTION(float extinction_ui, vec3 X) { return extinction_ui; }\n`; }
    if (shader.indexOf("VOLUME_EXTINCTION_MAX(")           == -1) { shader += `\n float VOLUME_EXTINCTION_MAX(float extinction_ui) { return extinction_ui; }\n`; }
    if (shader.indexOf("VOLUME_SCATTERING_COLOR(")         == -1) { shader += `\n vec3 VOLUME_SCATTERING_COLOR(vec3 scattering_color_ui, vec3 X) { return scattering_color_ui; }\n`; }
    if (shader.indexOf("VOLUME_ABSORPTION_COLOR(")         == -1) { shader += `\n vec3 VOLUME_ABSORPTION_COLOR(vec3 absorption_color_ui, vec3 X) { return absorption_color_ui; }\n`; }
    if (shader.indexOf("VOLUME_EMISSION(")                 == -1) { shader += `\n vec3 VOLUME_EMISSION(vec3 emission, vec3 X) { return emission; }\n`; }
    if (shader.indexOf("VOLUME_ANISOTROPY(")               == -1) { shader += `\n float VOLUME_ANISOTROPY(float anisotropy, vec3 X) { return anisotropy; }\n`; }

    let hasGeometry = (hasSurface || hasMetal || hasDielectric);
    if ( !(hasGeometry || hasVolume) )
    { 
        GLU.fail('Scene must define at least one of: SDF_SURFACE, SDF_METAL, SDF_DIELECTRIC, or SDF_VOLUME'); 
    }

    let hasNM = (hasSurfaceNM || hasMetalNM || hasDielectricNM);
    
    // Insert material GLSL code
    var dielectricObj = snelly.getLoadedDielectric(); if (dielectricObj == null) return;
    var metalObj      = snelly.getLoadedMetal();      if (metalObj == null) return;
    iorCodeDiele    = dielectricObj.ior();
    iorCodeMetal    = metalObj.ior();

    // Copy the current scene and material routines into the source code
    // of the trace fragment shader
    replacements = {};
    replacements.__SHADER__          = shader;
    replacements.__IOR_FUNC__        = iorCodeDiele + '\n' + iorCodeMetal;
    replacements.__MAX_MARCH_STEPS__ = Math.round(this.maxMarchSteps);
    replacements.__MAX_VOLUME_STEPS__ = Math.round(this.maxVolumeSteps);
    replacements.__MAX_BOUNCES__  = Math.round(this.maxBounces);
    replacements.__MAX_SCATTERS__ = Math.round(this.maxScatters);
    replacements.__MAX_SAMPLES_PER_FRAME__  = Math.round(this.maxSamplesPerFrame); 

    replacements.__DEFINES__ = '';
    if (hasSurface)    replacements.__DEFINES__ += '\n#define HAS_SURFACE\n';
    if (hasMetal)      replacements.__DEFINES__ += '\n#define HAS_METAL\n';
    if (hasDielectric) replacements.__DEFINES__ += '\n#define HAS_DIELECTRIC\n';
    if (hasVolume)     replacements.__DEFINES__ += '\n#define HAS_VOLUME\n';
    if (hasGeometry)   replacements.__DEFINES__ += '\n#define HAS_GEOMETRY\n';

    if (hasNM)           replacements.__DEFINES__ += '\n#define HAS_NORMALMAP\n';
    if (hasSurfaceNM)    replacements.__DEFINES__ += '\n#define HAS_SURFACE_NORMALMAP\n';
    if (hasMetalNM)      replacements.__DEFINES__ += '\n#define HAS_METAL_NORMALMAP\n';
    if (hasDielectricNM) replacements.__DEFINES__ += '\n#define HAS_DIELECTRIC_NORMALMAP\n';

    if (this.interactive) replacements.__DEFINES__ += '\n#define INTERACTIVE_MODE\n';

    // Compile pathtracer with different entry point according to mode.
    // Here shaderSources is a dict from name (e.g. "trace")
    // to a dict {v:vertexShaderSource, f:fragmentShaderSource}
    switch (this.renderMode)
    {
        case 'ao':
            this.aoProgram = new GLU.Shader('ao', this.shaderSources, replacements);
            break;
        case 'normals':
            this.normalsProgram = new GLU.Shader('normals', this.shaderSources, replacements);
            break;
        case 'firsthit':
            if (!hasSurface) 
            {
                GLU.fail('"firsthit" shader requires definition of SDF_SURFACE'); 
            }
            this.firsthitProgram = new GLU.Shader('firsthit', this.shaderSources, replacements);
            break;
        case 'pt':
        default:
            this.pathtraceAllProgram = new GLU.Shader('pathtracer', this.shaderSources, replacements);
            break;
    }

    this.pickProgram = new GLU.Shader('pick', this.shaderSources, replacements);

    // Tonemapping program
    this.tonemapProgram = new GLU.Shader('tonemapper', this.shaderSources, null);
}

Renderer.prototype.enabled = function()
{
    return this.enable;
}

Renderer.prototype.depthTestEnabled = function()
{
    return this.depthTest;
}

Renderer.prototype.updateSunDir = function()
{
    let latTheta = (90.0-this.sunLatitude) * Math.PI/180.0;
    let lonPhi = this.sunLongitude * Math.PI/180.0;
    let costheta = Math.cos(latTheta);
    let sintheta = Math.sin(latTheta);
    let cosphi = Math.cos(lonPhi);
    let sinphi = Math.sin(lonPhi);
    let x = sintheta * cosphi;
    let z = sintheta * sinphi;
    let y = costheta;
    this.sunDir = [x, y, z];
}

/**
* Adjust focal length to match the distance to the (first) object visible at a given location in the frame.
* @constructor
* @property {number} x - x coordinate of focus as a fraction of current frame width
* @property {number} y - y coordinate of focus as a fraction of current frame width
*/
Renderer.prototype.focusAt = function(xPick, yPick)
{
    let px = x * this._width;
    let py = y * this._height;
    this.pick(px, py);
}

Renderer.prototype.pick = function(xPick, yPick)
{
    if (!this.loaded) return;
    if (this.pickProgram==null) return;
    var sceneObj = snelly.getScene(); if (sceneObj==null) return;

    let gl = this.gl;
    gl.viewport(0, 0, 1, 1);

    let PROGRAM = this.pickProgram;
    PROGRAM.bind();

    // Upload current scene shader parameters
    if (typeof sceneObj.syncShader !== "undefined")
    {
        sceneObj.syncShader(snelly, this.pickProgram);
    }

    // sync camera info to shader
    var camera = snelly.getCamera();
    var camPos = camera.position.clone();
    var camDir = camera.getWorldDirection();
    var camUp = camera.up.clone();
    camUp.transformDirection( camera.matrixWorld );
    var camX = new THREE.Vector3();
    camX.crossVectors(camUp, camDir);
    PROGRAM.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
    PROGRAM.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
    PROGRAM.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
    PROGRAM.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
    PROGRAM.uniformF("camFovy", camera.fov);
    PROGRAM.uniformF("camAspect", camera.aspect);
    PROGRAM.uniformF("camAperture", snelly.lengthScale*Math.pow(1.4142135,camera.aperture));
    PROGRAM.uniformF("camFocalDistance", snelly.lengthScale*Math.pow(10.0,camera.focalDistance));
    PROGRAM.uniform2Fv("resolution", [this._width, this._height]);
    PROGRAM.uniformF("lengthScale", snelly.lengthScale);
    PROGRAM.uniformF("minLengthScale", snelly.minLengthScale);
    PROGRAM.uniformF("maxLengthScale", snelly.maxLengthScale);
    PROGRAM.uniformI("maxStepsIsMiss", Boolean(this.maxStepsIsMiss) ? 1 : 0);
    PROGRAM.uniform2Fv("mousePick", [xPick, yPick]); // picked pixel

    let fbo = new GLU.RenderTarget();
    fbo.bind();
    fbo.drawBuffers(1);
    gl.bindTexture(gl.TEXTURE_2D, null);

    let pickData = new Float32Array(4);
    this.pickTex = new GLU.Texture(1, 1, 4, true, false, true, pickData);
    fbo.attachTexture(this.pickTex, 0);

    // Trace pick ray
    this.quadVbo.bind();
    this.quadVbo.draw(PROGRAM, gl.TRIANGLE_FAN);

    // Read floating point hit distance output from pick fragment shader
    var pixels = new Float32Array(4);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, pixels);
    fbo.unbind();

    return {distance: pixels[0],
            material: pixels[1]};
}

Renderer.prototype.render = function()
{
    if (!this.loaded) return;
    var sceneObj = snelly.getScene(); if (sceneObj==null) return;
    if (snelly.getSpectra()==null) return;

    var timer_start = performance.now();

    let gl = this.gl;
    if (typeof sceneObj.preframeCallback != "undefined")
    {
        sceneObj.preframeCallback(snelly, gl);
    }

    gl.disable(gl.DEPTH_TEST);
    gl.viewport(0, 0, this._width, this._height);

    // Update expected sample count after this frame, based on current resolution and skip probability
    if (this.interactive)
    {
        this.numSamples += (1.0-this.skipProbability) * Math.round(this.maxSamplesPerFrame) * this._width * this._height;
    }
    else
    {   
        this.numSamples += Math.round(this.maxSamplesPerFrame) * this._width * this._height;
    }
    this.spp = this.numSamples / (this._width * this._height);

    // Update framebuffer only if we reached the requied SPP threshold for redraw
    if (this.spp >= this.minsSPPToRedraw)
    {
        gl.clearColor(0.0, 0.0, 0.0, 1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }

    ////////////////////////////////////////////////
    /// Pathtracing
    ////////////////////////////////////////////////

    let INTEGRATOR_PROGRAM = null;
    switch (this.renderMode)
    {
        case 'ao':       INTEGRATOR_PROGRAM = this.aoProgram;           break;
        case 'normals':  INTEGRATOR_PROGRAM = this.normalsProgram;      break;
        case 'firsthit': INTEGRATOR_PROGRAM = this.firsthitProgram;     break;
        case 'pt':
        default:         INTEGRATOR_PROGRAM = this.pathtraceAllProgram; break;
    }

    INTEGRATOR_PROGRAM.bind();

    // sync camera info to shader
    var camera = snelly.getCamera();
    var camPos = camera.position.clone();
    var camDir = camera.getWorldDirection();
    var camUp = camera.up.clone();
    camUp.transformDirection( camera.matrixWorld );
    var camX = new THREE.Vector3();
    camX.crossVectors(camUp, camDir);
    INTEGRATOR_PROGRAM.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
    INTEGRATOR_PROGRAM.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
    INTEGRATOR_PROGRAM.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
    INTEGRATOR_PROGRAM.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
    INTEGRATOR_PROGRAM.uniformF("camFovy", camera.fov);
    INTEGRATOR_PROGRAM.uniformF("camAspect", camera.aspect);
    INTEGRATOR_PROGRAM.uniformF("camAperture", snelly.lengthScale*Math.pow(1.4142135,camera.aperture));
    INTEGRATOR_PROGRAM.uniformF("camFocalDistance", snelly.lengthScale*Math.pow(10.0,camera.focalDistance));
    INTEGRATOR_PROGRAM.uniform2Fv("resolution", [this._width, this._height]);

    // Read wavelength -> XYZ table
    if (this.renderMode=='pt' || this.renderMode=='ao' || this.renderMode=='firsthit')
    {
        snelly.wavelengthToXYZ.bind(2);
        INTEGRATOR_PROGRAM.uniformTexture("WavelengthToXYZ", snelly.wavelengthToXYZ);
        snelly.emissionIcdf.bind(3);
        INTEGRATOR_PROGRAM.uniformTexture("ICDF", snelly.emissionIcdf);
    }

    // Pathtracing options
    INTEGRATOR_PROGRAM.uniformF("skyPower", this.skyPower);
    INTEGRATOR_PROGRAM.uniformF("sunPower", this.sunPower);
    INTEGRATOR_PROGRAM.uniformF("sunAngularSize", this.sunAngularSize);
    INTEGRATOR_PROGRAM.uniformF("sunLatitude", this.sunLatitude);
    INTEGRATOR_PROGRAM.uniformF("sunLongitude", this.sunLongitude);
    INTEGRATOR_PROGRAM.uniform3Fv("sunColor", this.sunColor);
    this.updateSunDir();
    INTEGRATOR_PROGRAM.uniform3Fv("sunDir", this.sunDir);
    INTEGRATOR_PROGRAM.uniformI("sunVisibleDirectly", this.sunVisibleDirectly);

    INTEGRATOR_PROGRAM.uniformF("radianceClamp", Math.pow(10.0, this.radianceClamp));
    INTEGRATOR_PROGRAM.uniformF("skipProbability", this.skipProbability);
    INTEGRATOR_PROGRAM.uniformF("lengthScale", Math.max(snelly.lengthScale, 1.0e-6));
    INTEGRATOR_PROGRAM.uniformF("maxLengthScale", Math.max(snelly.maxLengthScale, 1.0e-6));
    INTEGRATOR_PROGRAM.uniformF("minLengthScale", Math.max(snelly.minLengthScale, 1.0e-6));
    INTEGRATOR_PROGRAM.uniformF("shadowStrength", this.shadowStrength);
    INTEGRATOR_PROGRAM.uniformI("envMapVisible", Boolean(this.envMapVisible) ? 1 : 0);
    INTEGRATOR_PROGRAM.uniformF("envMapRotation", Math.min(Math.max(this.envMapRotation, 0.0), 360.0));
    INTEGRATOR_PROGRAM.uniformI("maxStepsIsMiss", Boolean(this.maxStepsIsMiss) ? 1 : 0);
    INTEGRATOR_PROGRAM.uniformI("wavelengthSamples", this.wavelengthSamples);
    
    // Attach radiance FBO
    this.fbo.bind();
    this.fbo.drawBuffers(2);
    var current = this.currentState;
    var next    = 1 - current;
    this.pathStates[current].bind(INTEGRATOR_PROGRAM); // Read data from the 'current' state
    this.pathStates[next].attach(this.fbo);            // Write data into the 'next' state

    // Bind env map if we have it
    INTEGRATOR_PROGRAM.uniformI("haveEnvMap", Boolean(this.envMap) ? 1 : 0);
    if (this.envMap != null)
    {
        gl.activeTexture(gl.TEXTURE0 + 6);
        gl.bindTexture(gl.TEXTURE_2D, this.envMap.tex);
        var id = gl.getUniformLocation(INTEGRATOR_PROGRAM.program, "envMap");
        gl.uniform1i(id, 6);
    }

    // Upload current scene shader parameters
    if (typeof sceneObj.syncShader !== "undefined")
    {
        sceneObj.syncShader(snelly, INTEGRATOR_PROGRAM);
    }

    // Upload material parameters
    snelly.materials.syncShader(INTEGRATOR_PROGRAM);

    // Trace one path per pixel
    gl.disable(gl.BLEND);
    this.quadVbo.bind();
    this.quadVbo.draw(INTEGRATOR_PROGRAM, gl.TRIANGLE_FAN);
    this.fbo.unbind();
    gl.bindTexture(gl.TEXTURE_2D, null);

    ////////////////////////////////////////////////
    /// Tonemapping / compositing
    ////////////////////////////////////////////////

    if (this.spp >= this.minsSPPToRedraw)
    {
        this.tonemapProgram.bind();
        var radianceTexCurrent = this.pathStates[next].radianceTex;
        radianceTexCurrent.bind(0);
        this.tonemapProgram.uniformTexture("Radiance", radianceTexCurrent);
        this.tonemapProgram.uniformF("exposure", this.exposure);
        this.tonemapProgram.uniformF("invGamma", 1.0/this.gamma);
        this.tonemapProgram.uniformF("contrast", this.contrast);
        this.tonemapProgram.uniformF("saturation", this.saturation);

        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
        gl.blendEquation(gl.FUNC_ADD);
        this.quadVbo.bind();
        this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);
    }

    // Ping-pong radiance buffers
    this.currentState = next;

    // Frame timing
    var timer_end = performance.now();
    var frame_time_ms = (timer_end - timer_start);
    if (this.numSamples==0)
    {
        this.frametime_measure_ms = frame_time_ms;
    }
    else
    {
        // apply exponential smoothing to frame time measurements
        var smoothing = 0.5;
        var prev_frame_time_ms = this.frametime_measure_ms;
        this.frametime_measure_ms = smoothing*frame_time_ms + (1.0-smoothing)*prev_frame_time_ms;
    }
    this.frametime_measure_ms = Math.max(1.0-6, this.frametime_measure_ms);

    // Update skip probability for next frame, to try to achieve interactive framerate
    let goalMs = Math.min(1.0e3, Math.max(1.0, 1.0e3/this.goalFPS));
    this.skipProbability = Math.max(0.0, 1.0-goalMs/this.frametime_measure_ms);

    this.numFramesSinceReset++;
    this.numFramesSinceInit++;

    if (typeof sceneObj.postframeCallback != "undefined")
    {
        sceneObj.postframeCallback(snelly, gl);
    }
}

Renderer.prototype.resize = function(width, height)
{
    this._width = width;
    this._height = height;

    this.fbo.unbind();
    // Two sets of radiance buffers, for read/write ping-pong
    this.pathStates = [new PathtracerState(this._width, this._height),
                       new PathtracerState(this._width, this._height)];

    this.quadVbo = this.createQuadVbo();
    this.fbo = new GLU.RenderTarget();
    this.reset(true);
}
