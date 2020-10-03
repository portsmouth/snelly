
/**
* Snelly is the global object providing access to all functionality in the system.
* @constructor
* @param {Scene} sceneObj - The user-defined scene
*/
var Snelly = function(sceneObj)
{
    console.log('[snelly] Snelly constructor');

    this.initialized = false;
    this.terminated = false;
    this.rendering = false;
    this.sceneObj = sceneObj;
    snelly = this;

    let container = document.getElementById("container");
    this.container = container;

    var render_canvas = document.getElementById('render-canvas');
    this.render_canvas = render_canvas;
    this.width = render_canvas.width;
    this.height = render_canvas.height;
    render_canvas.style.width = render_canvas.width;
    render_canvas.style.height = render_canvas.height;

    var text_canvas = document.getElementById('text-canvas');
    this.text_canvas = text_canvas;
    this.textCtx = text_canvas.getContext("2d");
    this.onSnellyLink = false;
    this.onUserLink = false;
    this.statusText = '';

    window.addEventListener( 'resize', this, false );

    // Setup THREE.js orbit camera
    var VIEW_ANGLE = 45; // @todo: fov should be under user control
    var ASPECT = this.width / this.height;
    var NEAR = 0.05;
    var FAR = 1000;
    this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
    this.camera.lookAt(new THREE.Vector3(0.0, 0.0, 0.0));
    this.camera.position.set(1.0, 1.0, 1.0);

    this.camControls = new THREE.OrbitControls(this.camera, this.container);
    this.camControls.zoomSpeed = 2.0;
    this.camControls.flySpeed = 0.01;
    this.camControls.addEventListener('change', camChanged);
    this.camControls.keyPanSpeed = 100.0;

    this.gui = null;
    this.guiVisible = true;

    // Instantiate materials
    this.materials = new Materials();

    // Instantiate distance field pathtracer
    this.pathtracer = new Renderer();
    this.auto_resize = true;

    // Spectrum initialization
    this.spectra = {}
    this.SPECTRUM_SAMPLES = 1024;
    this.spectrumObj = null;
    this.LAMBDA_MIN = 390.0;
    this.LAMBDA_MAX = 750.0;
    var wToXYZ = wavelengthToXYZTable();
    this.wavelengthToXYZ = new GLU.Texture(wToXYZ.length/4, 1, 4, true, true, true, wToXYZ);
    this.emissionIcdf    = new GLU.Texture(4*this.SPECTRUM_SAMPLES, 1, 1, true, true, true, null);

    // Allow user to programmatically initialize the camera, materials, and renderer
    this.initScene();

    // Do initial resize:
    this.resize();

    // Create dat gui
    this.gui = new GUI(this.guiVisible);

    // Setup keypress and mouse events
    window.addEventListener( 'mousemove', this, false );
    window.addEventListener( 'mousedown', this, false );
    window.addEventListener( 'mouseup',   this, false );
    window.addEventListener( 'contextmenu',   this, false );
    window.addEventListener( 'click', this, false );
    window.addEventListener( 'keydown', this, false );

    this.reset_time = performance.now();
    this.initialized = true;
}

/**
* Returns the current version number of the snelly system, in the format [1, 2, 3] (i.e. major, minor, patch version)
*  @returns {Array}
*/
Snelly.prototype.getVersion = function()
{
    return [1, 11, 0];
}

Snelly.prototype.handleEvent = function(event)
{
    switch (event.type)
    {
        case 'resize':      this.resize();  break;
        case 'mousemove':   this.onDocumentMouseMove(event);  break;
        case 'mousedown':   this.onDocumentMouseDown(event);  break;
        case 'mouseup':     this.onDocumentMouseUp(event);    break;
        case 'contextmenu': this.onDocumentRightClick(event); break;
        case 'click':       this.onClick(event);  break;
        case 'keydown':     this.onkeydown(event);  break;
    }
}

/**
* Access to the Renderer object
*  @returns {Renderer}
*/
Snelly.prototype.getRenderer = function()
{
    return this.pathtracer;
}

/**
* Access to the GUI object
*  @returns {GUI}
*/
Snelly.prototype.getGUI = function()
{
    return this.gui;
}

/**
* Access to the camera object
* @returns {THREE.PerspectiveCamera}.
*/
Snelly.prototype.getCamera = function()
{
    return this.camera;
}

/**
* Access to the camera controller object
* @returns {THREE.OrbitControls}
*/
Snelly.prototype.getControls = function()
{
    return this.camControls;
}

/**
* Programmatically show or hide the dat.GUI UI
* @param {Boolean} show - toggle
*/
Snelly.prototype.showGUI = function(show)
{
    this.guiVisible = show;
}

/**
* Specify arbitrary status text (one line only) to display in the lower right of the viewport
* @param {Boolean} statusText - text to display
*/
Snelly.prototype.setStatus = function(statusText)
{
    this.statusText = statusText;
}

//
// Scene management
//

Snelly.prototype.getScene = function()
{
    return this.sceneObj;
}

/**
 * @returns {WebGLRenderingContext} The webGL context
 */
Snelly.prototype.getGLContext = function()
{
    return GLU.gl;
}


Snelly.prototype.initScene = function()
{
    console.warn('[snelly] Snelly.prototype.initScene');

    if (typeof this.sceneObj.shader == "undefined")
    {
        GLU.fail('Scene must define a "shader" function!');
    }

    // Obtain length scales, if specified
    this.lengthScale = 1.0;
    if (typeof this.sceneObj.getLengthScale !== "undefined")
    {
        this.lengthScale = this.sceneObj.getLengthScale();
    }

    this.minLengthScale = 1.0e-4;
    if (typeof this.sceneObj.getMinLengthScale !== "undefined")
    {
        this.minLengthScale = this.sceneObj.getMinLengthScale();
    }
    this.minLengthScale = Math.max(1.0e-12, this.minLengthScale);  // sanity

    this.maxLengthScale = 1.0e2;
    if (typeof this.sceneObj.getMaxLengthScale !== "undefined")
    {
        this.maxLengthScale = this.sceneObj.getMaxLengthScale();
    }
    this.maxLengthScale = Math.min(1.0e9, this.maxLengthScale); // sanity

    // Set initial default camera position and target based on max scale
    let po = 2.0*this.lengthScale;
    this.camera.position.set(po, po, po);
    this.camControls.target.set(0.0, 0.0, 0.0);

    this.camera.aperture      = -10.0; // logarithmic (base sqrt(2)), relative to length scale
    this.camera.focalDistance =  1.0; // logarithmic (base 10), relative to length scale

    // Call user-defined init function
    if (typeof this.sceneObj.init !== "undefined")
    {
        this.sceneObj.init(this);
    }

    // Read optional scene name and URL, if provided
    this.sceneName = '';
    if (typeof this.sceneObj.getName !== "undefined")
    {
        this.sceneName = this.sceneObj.getName();
    }

    this.sceneURL = '';
    if (typeof this.sceneObj.getURL !== "undefined")
    {
        this.sceneURL = this.sceneObj.getURL();
    }

    // cache initial camera position to allow reset on 'F'
    this.initial_camera_position = new THREE.Vector3();
    this.initial_camera_position.copy(this.camera.position);
    this.initial_camera_target = new THREE.Vector3();
    this.initial_camera_target.copy(this.camControls.target);

    // Compile GLSL shaders
    this.pathtracer.compileShaders();

    // Fix renderer to width & height, if they were specified
    if ((typeof this.pathtracer.width!=="undefined") && (typeof this.pathtracer.height!=="undefined"))
    {
        this.auto_resize = false;
        this._resize(this.pathtracer.width, this.pathtracer.height);
    }

    // Set up blackbody spectrum sampling
    this.spectra = [];
    this.addSpectrum( new BlackbodySpectrum("blackbody", "Blackbody spectrum", this.pathtracer.skyTemperature) );
    this.loadSpectrum("blackbody");
    //this.addSpectrum( new FlatSpectrum("flat", "Flat spectrum", 390.0, 750.0) );
    //this.loadSpectrum("flat");

    // Camera setup
    this.camControls.update();
    this.reset(false);
}


Snelly.prototype.dumpScene = function()
{
    console.warn('[snelly] Snelly.prototype.dumpScene');

    let camera = this.camera;
    let controls = this.camControls;
    let renderer = this.pathtracer;
    let materials = this.materials;

    var code = `
{
    /******* copy-pasted console output on 'O', begin *******/\n`;
    code += `
    let renderer  = snelly.getRenderer();
    let camera    = snelly.getCamera();
    let controls  = snelly.getControls();
    let materials = snelly.getMaterials();
    `;

    if (typeof this.sceneObj.initGenerator !== "undefined")
    {
        code += this.sceneObj.initGenerator();
    }

    code += this.guiVisible ? `\nsnelly.showGUI(true);\n` : `\nsnelly.showGUI(false);\n`;

    code += `
    /** Camera settings **/
    camera.fov = ${camera.fov};
    camera.aperture = ${camera.aperture};
    camera.focalDistance = ${camera.focalDistance};
    camera.up.set(${camera.up.x}, ${camera.up.y}, ${camera.up.z});
    camera.position.set(${camera.position.x}, ${camera.position.y}, ${camera.position.z});
    controls.target.set(${controls.target.x}, ${controls.target.y}, ${controls.target.z});

    /** Renderer settings **/
    // General rendering settings
    renderer.renderMode = '${renderer.renderMode}';
    renderer.dispersive = ${renderer.dispersive};
    renderer.maxSamplesPerFrame = ${renderer.maxSamplesPerFrame};
    renderer.maxSpp = ${renderer.maxSpp};
    renderer.maxBounces = ${renderer.maxBounces};
    renderer.maxAtmosphereScatters = ${renderer.maxAtmosphereScatters};
    renderer.maxMarchSteps = ${renderer.maxMarchSteps};
    renderer.maxStepsIsMiss = ${renderer.maxStepsIsMiss};
    renderer.interactive = ${renderer.interactive};
    renderer.goalFPS = ${renderer.goalFPS};
    renderer.minsSPPToRedraw = ${renderer.minsSPPToRedraw};
    renderer.filterRadius = ${renderer.filterRadius};
    renderer.radianceClamp = ${renderer.radianceClamp};
    renderer.wavelengthSamples = ${renderer.wavelengthSamples};
    renderer.shadowStrength = ${renderer.shadowStrength};
    // Tone-mapping
    renderer.exposure = ${renderer.exposure};
    renderer.gamma = ${renderer.gamma};
    renderer.contrast = ${renderer.contrast};
    renderer.saturation = ${renderer.saturation};
    // Lights
        // sky light
        renderer.skyPower = ${renderer.skyPower};
        renderer.skyTintUp = [${renderer.skyTintUp[0]}, ${renderer.skyTintUp[1]}, ${renderer.skyTintUp[2]}];
        renderer.skyTintDown = [${renderer.skyTintDown[0]}, ${renderer.skyTintDown[1]}, ${renderer.skyTintDown[2]}];
        renderer.envMapVisible = ${renderer.envMapVisible};
        renderer.envMapPhiRotation = ${renderer.envMapPhiRotation};
        renderer.envMapThetaRotation = ${renderer.envMapThetaRotation};
        renderer.envMapTransitionAngle = ${renderer.envMapTransitionAngle};
        // sun light
        renderer.sunPower = ${renderer.sunPower};
        renderer.sunColor = [${renderer.sunColor[0]}, ${renderer.sunColor[1]}, ${renderer.sunColor[2]}];
        renderer.sunAngularSize = ${renderer.sunAngularSize};
        renderer.sunLatitude = ${renderer.sunLatitude};
        renderer.sunLongitude = ${renderer.sunLongitude};
        renderer.sunVisibleDirectly = ${renderer.sunVisibleDirectly};
        // sphere light
        renderer.sphereLightPosition = [${renderer.sphereLightPosition[0]}, ${renderer.sphereLightPosition[1]}, ${renderer.sphereLightPosition[2]}];
        renderer.sphereLightRadius = ${renderer.sphereLightRadius};
        renderer.sphereLightPower = ${renderer.sphereLightPower};
        renderer.sphereLightColor = [${renderer.sphereLightColor[0]}, ${renderer.sphereLightColor[1]}, ${renderer.sphereLightColor[2]}];
`;

    code += `
    /** Material settings **/`;
    var shader = this.sceneObj.shader();
    if (shader.indexOf("SDF_SURFACE(") != -1)
    {
        code += `
    let surface = materials.loadSurface();`;
        code += materials.loadSurface().repr();
    }
    if (shader.indexOf("SDF_METAL(") != -1)
    {
        code += `
    let metal = materials.loadMetal('${materials.getLoadedMetal().getName()}');`;
        code += materials.loadMetal(materials.getLoadedMetal().getName()).repr();
    }
    if (shader.indexOf("SDF_DIELECTRIC(") != -1)
    {
        code += `
    let dielectric = materials.loadDielectric('${materials.getLoadedDielectric().getName()}');`;
        code += materials.loadDielectric(materials.getLoadedDielectric().getName()).repr();
    }
    code += `
    let volume = materials.loadVolume();`;
    code += materials.loadVolume().repr();

    code += `
    /******* copy-pasted console output on 'O', end *******/
}`;
    return code;
}

// emission spectrum management
Snelly.prototype.addSpectrum = function(spectrumObj)
{
    this.spectra[spectrumObj.getName()] = spectrumObj;
}

Snelly.prototype.getSpectra = function()
{
    return this.spectra;
}

Snelly.prototype.loadSpectrum = function(spectrumName)
{
    this.spectrumObj = this.spectra[spectrumName];
    var inverseCDF = this.spectrumObj.inverseCDF(this.LAMBDA_MIN, this.LAMBDA_MAX, this.SPECTRUM_SAMPLES);
    this.emissionIcdf.bind(0);
    this.emissionIcdf.copy(inverseCDF);
    this.reset(true);
}

Snelly.prototype.getLoadedSpectrum = function()
{
    return this.spectrumObj;
}


// Material management

/**
* Get materials object
 * @returns {Materials}
*/
Snelly.prototype.getMaterials = function()
{
    return this.materials;
}

Snelly.prototype.getDielectrics = function()
{
    return this.materials.getDielectrics();
}

Snelly.prototype.getMetals = function()
{
    return this.materials.getMetals();
}

Snelly.prototype.loadDielectric = function(dielectricName)
{
    this.materials.loadDielectric(dielectricName);
    this.reset();
}

Snelly.prototype.loadMetal = function(metalName)
{
    this.materials.loadMetal(metalName);
    this.reset();
}

Snelly.prototype.getLoadedDielectric = function()
{
    return this.materials.getLoadedDielectric();
}

Snelly.prototype.getLoadedMetal = function()
{
    return this.materials.getLoadedMetal();
}

/**
* Get Surface object
 * @returns {Surface}
*/
Snelly.prototype.getSurface = function()
{
    return this.materials.loadSurface();
}

/**
* Get Volume object
 * @returns {Volume}
*/
Snelly.prototype.getVolume = function()
{
    return this.materials.loadVolume();
}

// Renderer reset on camera or other parameters update
Snelly.prototype.reset = function(no_recompile = false)
{
    if (!this.initialized || this.terminated) return;
    this.reset_time = performance.now();
    this.pathtracer.reset(no_recompile);
}
   
// Render all
Snelly.prototype.render = function()
{
    if (!this.initialized || this.terminated) return;
    if (this.sceneObj == null) return;
    this.rendering = true;

    let isReady = true;
    if (typeof this.sceneObj.isReady !== "undefined")
    {
        isReady = this.sceneObj.isReady(this);
    }

    if (isReady)
    {
        // Render distance field surface via pathtracing
        this.pathtracer.render();
    }
    else
    {
        const gl = GLU.gl;
        gl.clearColor(0.0, 0.0, 0.0, 1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }

    // Update HUD text canvas
    this.textCtx.textAlign = "left";   	// This determines the alignment of text, e.g. left, center, right
    this.textCtx.textBaseline = "middle";	// This determines the baseline of the text, e.g. top, middle, bottom
    this.textCtx.font = '12px monospace';	// This determines the size of the text and the font family used
    this.textCtx.clearRect(0, 0, this.textCtx.canvas.width, this.textCtx.canvas.height);
    this.textCtx.globalAlpha = 0.95;
    this.textCtx.strokeStyle = 'black';
    this.textCtx.lineWidth  = 2;
    if (this.guiVisible)
    {
          if (this.onSnellyLink) this.textCtx.fillStyle = "#ff5500";
          else                   this.textCtx.fillStyle = "#ffff00";
          let ver = this.getVersion();
          this.textCtx.strokeText('Snelly renderer v'+ver[0]+'.'+ver[1]+'.'+ver[2], 14, 20);
          this.textCtx.fillText('Snelly renderer v'+ver[0]+'.'+ver[1]+'.'+ver[2], 14, 20);
          
          this.textCtx.fillStyle = "#ffccaa";
          this.textCtx.strokeText('spp: ' + (this.pathtracer.spp).toPrecision(3), 14, 35);
          this.textCtx.fillText('spp: ' + (this.pathtracer.spp).toPrecision(3), 14, 35);
          if (this.sceneName != '')
          {
              this.textCtx.fillStyle = "#ffaa22";
              this.textCtx.strokeText(this.sceneName, 14, this.textCtx.canvas.height-25);
              this.textCtx.fillText(this.sceneName, 14, this.textCtx.canvas.height-25);
          }
        if (this.sceneURL != '')
        {
            if (this.onUserLink) this.textCtx.fillStyle = "#aaccff";
            else               this.textCtx.fillStyle = "#55aaff";
            this.textCtx.strokeText(this.sceneURL, 14, this.textCtx.canvas.height-40);
            this.textCtx.fillText(this.sceneURL, 14, this.textCtx.canvas.height-40);
        }
        if (this.statusText != '')
        {
            let since_reset_time_ms = performance.now() - this.reset_time;
            let phase = Math.abs(Math.cos(2.0*since_reset_time_ms/1000.0));
            let r = Math.round(32+(255-32)*phase);
            let g = Math.round(32+(128-32)*phase);
            this.textCtx.fillStyle = "rgb(" + r + ", " + g + ", 0)"; 
            let x = this.textCtx.canvas.width - 12*this.statusText.length;
            this.textCtx.strokeText(this.statusText, x, this.textCtx.canvas.height-25);
            this.textCtx.fillText(this.statusText, x, this.textCtx.canvas.height-25);
        }
    }

    this.rendering = false;
}

Snelly.prototype._resize = function(width, height)
{
    this.width = width;
    this.height = height;

    let render_canvas = this.render_canvas;
    render_canvas.width  = width;
    render_canvas.height = height;
    render_canvas.style.width = width;
    render_canvas.style.height = height;

    var text_canvas = this.text_canvas;
    text_canvas.width  = width;
    text_canvas.height = height

    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.camControls.update();

    this.pathtracer.resize(width, height);
}

Snelly.prototype.resize = function()
{
    console.warn('[snelly] Snelly.prototype.resize');

    if (this.terminated) return;
    if (this.auto_resize)
    {
        // If no explicit renderer size was set by user, resizing the browser window
        // resizes the render itself to match.
        let width = window.innerWidth;
        let height = window.innerHeight;
        this._resize(width, height);
        if (this.initialized)
            this.render();
    }
    else
    {
        // Otherwise if the user set a fixed renderer resolution, we scale the resultant render
        // to fit into the current window with preserved aspect ratio:
        let render_canvas = this.render_canvas;
        let window_width = window.innerWidth;
        let window_height = window.innerHeight;
        let render_aspect = render_canvas.width / render_canvas.height;
        let window_aspect = window_width / window_height;
        if (render_aspect > window_aspect)
        {
            render_canvas.style.width = window_width;
            render_canvas.style.height = window_width / render_aspect;
        }
        else
        {
            render_canvas.style.width = window_height * render_aspect;
            render_canvas.style.height = window_height;
        }
        var text_canvas = this.text_canvas;
        text_canvas.width = window_width;
        text_canvas.height = window_height;
    }
}


/**
*
* @returns {number} - the minimum texture unit for user supplied textures in the shader
*/
Snelly.prototype.getUserTextureUnitStart = function()
{
    return 7;
}

Snelly.prototype.onClick = function(event)
{
    if (this.onSnellyLink)
    {
        window.open("https://github.com/portsmouth/snelly");
    }
    if (this.onUserLink)
    {
        window.open(this.sceneURL);
    }
    event.preventDefault();
}

Snelly.prototype.onDocumentMouseMove = function(event)
{
    // Check whether user is trying to click the Snelly home link, or user link
    var textCtx = this.textCtx;
    var x = event.pageX;
    var y = event.pageY;
    let linkWidth = this.textCtx.measureText('Snelly renderer vX.X.X').width;
    if (x>14 && x<14+linkWidth && y>15 && y<25) this.onSnellyLink = true;
    else this.onSnellyLink = false;
    if (this.sceneURL != '')
    {
        linkWidth = this.textCtx.measureText(this.sceneURL).width;
        if (x>14 && x<14+linkWidth && y>this.height-45 && y<this.height-35) this.onUserLink = true;
        else this.onUserLink = false;
    }

    this.camControls.update();
}

Snelly.prototype.onDocumentMouseDown = function(event)
{
    this.camControls.update();
}

Snelly.prototype.onDocumentMouseUp = function(event)
{
    this.camControls.update();
}

Snelly.prototype.onDocumentRightClick = function(event)
{
    if (event.altKey) return; // don't pick if alt-right-clicking (panning)
    let render_canvas = this.render_canvas;

    // map pixel picked on window to pixel of render buffer
    let rsw = parseInt(render_canvas.style.width, 10);
    let rsh = parseInt(render_canvas.style.height, 10);
    let wxPick = Math.min(event.clientX, rsw);
    let wyPick = Math.max(0, rsh - event.clientY);

    let rxPick = (wxPick/rsw) * render_canvas.width;
    let ryPick = (wyPick/rsh) * render_canvas.height;
    var pickData = this.pathtracer.pick(rxPick, ryPick);
    let fd = Math.max(1.0e-6, Math.abs(pickData.distance));
    this.camera.focalDistance = Math.log10(fd/this.lengthScale);
    this.gui.sync();
    this.reset(true);
}

Snelly.prototype.onkeydown = function(event)
{
    var charCode = (event.which) ? event.which : event.keyCode;
    switch (charCode)
    {
        case 122: // F11 key: go fullscreen
            var element	= document.body;
            if      ( 'webkitCancelFullScreen' in document ) element.webkitRequestFullScreen();
            else if ( 'mozCancelFullScreen'    in document ) element.mozRequestFullScreen();
            else console.assert(false);
            break;

        case 70: // F key: reset cam  (@todo)
            this.camera.position.copy(this.initial_camera_position);
            this.camControls.target.copy(this.initial_camera_target);
            this.reset(true);
            break;

        case 82: // R key: reset scene
            this.initScene();
            break;

        case 72: // H key: toggle hide/show dat gui
            this.guiVisible = !this.guiVisible;
            snelly.getGUI().toggleHide();
            break;

        case 79: // O key: output scene settings code to console
            let code = this.dumpScene();
            console.log(code);
            break;

        case 80: // P key: save current image to disk
        {
            var currentdate = new Date(); 
            var datetime = currentdate.getDate() + "-" + (currentdate.getMonth()+1)  + "-" + currentdate.getFullYear() + "_"  
                         + currentdate.getHours() + "-" + currentdate.getMinutes() + "-" + currentdate.getSeconds();
            let filename = `snelly-screenshot-${datetime}.png`;
            let link = document.createElement('a');
            link.download = filename;
            this.render_canvas.toBlob(function(blob){
                    link.href = URL.createObjectURL(blob);
                    var event = new MouseEvent('click');
                    link.dispatchEvent(event);
                },'image/png', 1);
            break;
        }

        case 87: // W key: cam forward
        {
            let toTarget = new THREE.Vector3();
            toTarget.copy(this.camControls.target);
            toTarget.sub(this.camera.position);
            let distToTarget = toTarget.length();
            toTarget.normalize();
            var move = new THREE.Vector3();
            move.copy(toTarget);
            move.multiplyScalar(this.camControls.flySpeed*distToTarget);
            this.camera.position.add(move);
            this.camControls.target.add(move);
            this.reset(true);
            break;
        }

        case 65: // A key: cam left
        {
            let toTarget = new THREE.Vector3();
            toTarget.copy(this.camControls.target);
            toTarget.sub(this.camera.position);
            let distToTarget = toTarget.length();
            var localX = new THREE.Vector3(1.0, 0.0, 0.0);
            var worldX = localX.transformDirection( this.camera.matrix );
            var move = new THREE.Vector3();
            move.copy(worldX);
            move.multiplyScalar(-this.camControls.flySpeed*distToTarget);
            this.camera.position.add(move);
            this.camControls.target.add(move);
            this.reset(true);
            break;
        }

        case 83: // S key: cam back
        {
            let toTarget = new THREE.Vector3();
            toTarget.copy(this.camControls.target);
            toTarget.sub(this.camera.position);
            let distToTarget = toTarget.length();
            toTarget.normalize();
            var move = new THREE.Vector3();
            move.copy(toTarget);
            move.multiplyScalar(-this.camControls.flySpeed*distToTarget);
            this.camera.position.add(move);
            this.camControls.target.add(move);
            this.reset(true);
            break;
        }

        case 68: // S key: cam right
        {
            let toTarget = new THREE.Vector3();
            toTarget.copy(this.camControls.target);
            toTarget.sub(this.camera.position);
            let distToTarget = toTarget.length();
            var localX = new THREE.Vector3(1.0, 0.0, 0.0);
            var worldX = localX.transformDirection( this.camera.matrix );
            var move = new THREE.Vector3();
            move.copy(worldX);
            move.multiplyScalar(this.camControls.flySpeed*distToTarget);
            this.camera.position.add(move);
            this.camControls.target.add(move);
            this.reset(true);
            break;
        }
    }

    // Call user keydown callback
    if (typeof this.sceneObj.onkeydownCallback !== "undefined")
    {
        this.sceneObj.onkeydownCallback(event, snelly, GLU.gl);
    }
}

function camChanged()
{
    if (!snelly.rendering)
    {
        var no_recompile = true;
        snelly.reset(no_recompile);
    }
}
