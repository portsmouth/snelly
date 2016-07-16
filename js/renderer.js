


//////////////////////////////////////////////////////////////////////
// Laser pointer UI
//////////////////////////////////////////////////////////////////////

var LaserPointer = function(glRenderer, glScene, glCamera) 
{
	// Main group of all objects moving rigidly with the laser pointer:
	this.group = new THREE.Object3D();
	var group = this.group;

	this.camera = glCamera;
	this.raycaster = new THREE.Raycaster();

	this.mouse     = new THREE.Vector2();
	this.mousePrev = new THREE.Vector2();
	this.mousedownOffset = new THREE.Vector3();

	this.INTERSECTED = null;
	this.SELECTED = null;

	var housingGeo      = new THREE.CylinderGeometry( 0.2, 0.2, 2, 32 );
	var housingMaterial = new THREE.MeshPhongMaterial( { color: 0xdddddd, specular: 0x999999, shininess: 60 } );
	var housingObj      = new THREE.Mesh( housingGeo, housingMaterial  );
	group.add(housingObj);

	///////////////////////////////////////////////////////////
	// Intersection geometry for handle controls (invisible)
	///////////////////////////////////////////////////////////

	this.intersectionHandleGroup = new THREE.Object3D();
	var intersectionHandleGroup = this.intersectionHandleGroup;

	var rI = 0.2;
	var hI = 8.0;
	var RI = 2.0;

	// x-handle
	var xHandleIntersectionGeo       = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var xHandleIntersectionMaterial  = new THREE.MeshPhongMaterial( { color: 'red', specular: 0x999999, shininess: 60, visible: true } );
	var xHandleIntersectionObj       = new THREE.Mesh( xHandleIntersectionGeo, xHandleIntersectionMaterial  );
	xHandleIntersectionObj.rotation.z = 0.5*Math.PI;
	xHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(xHandleIntersectionObj);

	// y-handle
	var yHandleIntersectionGeo        = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var yHandleIntersectionMaterial   = new THREE.MeshPhongMaterial( { color: 'green', specular: 0x999999, shininess: 60, visible: true } );
	var yHandleIntersectionObj        = new THREE.Mesh( yHandleIntersectionGeo, yHandleIntersectionMaterial  );
	yHandleIntersectionObj.rotation.y = 0.5*Math.PI;
	yHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(yHandleIntersectionObj);

	// z-handle
	var zHandleIntersectionGeo        = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var zHandleIntersectionMaterial   = new THREE.MeshPhongMaterial( { color: 'blue', specular: 0x999999, shininess: 30, visible: true } );
	var zHandleIntersectionObj        = new THREE.Mesh( zHandleIntersectionGeo, zHandleIntersectionMaterial  );
	zHandleIntersectionObj.rotation.x = 0.5*Math.PI;
	zHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(zHandleIntersectionObj);

	// x-rot-handle
	var xRotHandleIntersectionGeo      = new THREE.TorusGeometry( RI, rI, 32, 100 );
	var xRotHandleIntersectionMaterial = new THREE.MeshPhongMaterial( { color: 'red', specular: 0x999999, shininess: 30, visible: true } );
	var xRotHandleIntersectionObj      = new THREE.Mesh( xRotHandleIntersectionGeo, xRotHandleIntersectionMaterial  );
	xRotHandleIntersectionObj.rotation.y = 0.5*Math.PI;
	intersectionHandleGroup.add(xRotHandleIntersectionObj);

	// z-rot-handle
	var zRotHandleIntersectionGeo      = new THREE.TorusGeometry( RI, rI, 32, 100 );
	var zRotHandleIntersectionMaterial = new THREE.MeshPhongMaterial( { color: 'blue', specular: 0x999999, shininess: 30, visible: true } );
	var zRotHandleIntersectionObj      = new THREE.Mesh( zRotHandleIntersectionGeo, zRotHandleIntersectionMaterial  );
	zRotHandleIntersectionObj.rotation.y = 0.0;
	intersectionHandleGroup.add(zRotHandleIntersectionObj);

	group.add(intersectionHandleGroup);


	//////////////////////////////////////////////////////////////////////////////
	// Rendered geometry for handle controls (thin, visually appealing lines)
	//////////////////////////////////////////////////////////////////////////////
	
	this.renderHandleGroup = new THREE.Object3D();
	var renderHandleGroup = this.renderHandleGroup;

	var rR = 0.05;
	var hR = 8.0;
	var RR = 2.0;

	// x-handle
	var xHandleRenderGeo      = new THREE.CylinderGeometry( rR, rR, hR, 32 );
	var xHandleRenderMaterial = new THREE.MeshPhongMaterial( { color: 'red', specular: 0x999999, shininess: 60 } );
	var xHandleRenderObj      = new THREE.Mesh( xHandleRenderGeo, xHandleRenderMaterial  );
	xHandleRenderObj.rotation.z = 0.5*Math.PI;
	xHandleRenderObj.position.x = -0.0;
	renderHandleGroup.add(xHandleRenderObj);

	// y-handle
	var yHandleRenderGeo      = new THREE.CylinderGeometry( rR, rR, hR, 32 );
	var yHandleRenderMaterial = new THREE.MeshPhongMaterial( { color: 'green', specular: 0x999999, shininess: 60 } );
	var yHandleRenderObj      = new THREE.Mesh( yHandleRenderGeo, yHandleRenderMaterial  );
	yHandleRenderObj.rotation.y = 0.5*Math.PI;
	yHandleRenderObj.position.x = -0.0;
	renderHandleGroup.add(yHandleRenderObj);

	// z-handle
	var zHandleRenderGeo      = new THREE.CylinderGeometry( rR, rR, hR, 32 );
	var zHandleRenderMaterial = new THREE.MeshPhongMaterial( { color: 'blue', specular: 0x999999, shininess: 30 } );
	var zHandleRenderObj      = new THREE.Mesh( zHandleRenderGeo, zHandleRenderMaterial  );
	zHandleRenderObj.rotation.x = 0.5*Math.PI;
	zHandleRenderObj.position.x = -0.0;
	renderHandleGroup.add(zHandleRenderObj);

	group.add(renderHandleGroup);

	// x-rot-handle
	var xRotHandleRenderGeo      = new THREE.TorusGeometry( RR, rR, 32, 100 );
	var xRotHandleMaterial = new THREE.MeshPhongMaterial( { color: 'red', specular: 0x999999, shininess: 30 } );
	var xRotHandleRenderObj = new THREE.Mesh( xRotHandleRenderGeo, xRotHandleMaterial  );
	xRotHandleRenderObj.rotation.y = 0.5*Math.PI;
	renderHandleGroup.add(xRotHandleRenderObj);

	// z-rot-handle
	var zRotHandleRenderGeo      = new THREE.TorusGeometry( RR, rR, 32, 100 );
	var zRotHandleMaterial = new THREE.MeshPhongMaterial( { color: 'blue', specular: 0x999999, shininess: 30 } );
	var zRotHandleRenderObj = new THREE.Mesh( zRotHandleRenderGeo, zRotHandleMaterial  );
	zRotHandleRenderObj.rotation.y = 0.0;
	renderHandleGroup.add(zRotHandleRenderObj);

	///////////////////////////////////////////////////////////
	// proxy to represent the emission point of the laser
	///////////////////////////////////////////////////////////

	var proxyGeo      = new THREE.SphereGeometry(0.1);
	var proxyMaterial = new THREE.MeshPhongMaterial( { color: 0xff0000, specular: 0x999900, shininess: 10 } );
	var proxyObj      = new THREE.Mesh( proxyGeo, proxyMaterial  );
	proxyObj.position.y = 1.0;
	group.add(proxyObj);
	group.position.x = 0.0;

	///////////////////////////////////////////////////////////
	// back plane to implement dragging 
	///////////////////////////////////////////////////////////

	var backPlaneObj =  new THREE.Mesh(
						new THREE.PlaneGeometry(2000, 2000),
						new THREE.MeshBasicMaterial( { visible: false } )
				);
	backPlaneObj.position.x = 0.0;
	backPlaneObj.position.y = 0.0;

	this.X = new THREE.Vector3(1, 0, 0);
	this.Y = new THREE.Vector3(0, 1, 0);
	this.Z = new THREE.Vector3(0, 0, 1);

	this.objects = {

		"housing": housingObj,
		"proxy": proxyObj,

		"xHandleIntersection": xHandleIntersectionObj,
		"yHandleIntersection": yHandleIntersectionObj,
		"zHandleIntersection": zHandleIntersectionObj,
		"xRotHandleIntersection": xRotHandleIntersectionObj,
		"zRotHandleIntersection": zRotHandleIntersectionObj,
		"intersectionHandleGroup" : intersectionHandleGroup,

		"xHandleRender": xHandleRenderObj,
		"yHandleRender": yHandleRenderObj,
		"zHandleRender": zHandleRenderObj,
		"xRotHandleRender": xRotHandleRenderObj,
		"zRotHandleRender": zRotHandleRenderObj,
		"renderHandleGroup" : renderHandleGroup,
		
		"backPlane": backPlaneObj,
		"group": group
	};

	// Set these names on the objects
	for (var key in this.objects) 
	{
		if (this.objects.hasOwnProperty(key)) 
		{
			this.objects[key].name = key;
		}
	}

	var intersectionHandleNames = [ "xHandleIntersection", 
									"yHandleIntersection", 
									"zHandleIntersection",
									"xRotHandleIntersection", 
									"zRotHandleIntersection" ];
	this.intersectionHandleObjects = [];
	for (var n=0; n<intersectionHandleNames.length; n++) 
	{
		this.intersectionHandleObjects.push(this.objects[intersectionHandleNames[n]]);
	}

	var renderHandleNames = [ "xHandleRender", 
							  "yHandleRender", 
							  "zHandleRender",
							  "xRotHandleRender", 
							  "zRotHandleRender" ];
	this.renderHandleObjects = [];
	for (var n=0; n<renderHandleNames.length; n++) 
	{
		this.renderHandleObjects.push(this.objects[renderHandleNames[n]]);
	}

	this.intersectionHandleNameToRenderHandleName = { "xHandleIntersection": "xHandleRender",
													  "yHandleIntersection": "yHandleRender",
													  "zHandleIntersection": "zHandleRender",
													  "xRotHandleIntersection": "xRotHandleRender",
													  "zRotHandleIntersection": "zRotHandleRender" };
	this.glScene = glScene;
	this.glRenderer = glRenderer;
	this.glCamera = glCamera;

	// Finally, add the laser pointer objects (and backplane for raycasting) to the scene
	glScene.add(group);
	glScene.add(backPlaneObj);
}


LaserPointer.prototype.render = function()
{
	// Ensure laser handles are scaled according to current camera position
	// (to be roughly constant size in screen space)
	var group = this.objects["group"];

	var camDist = new THREE.Vector3();
	camDist.copy(this.objects["group"].position).sub(this.camera.position);
	var C = 0.05*camDist.length();

	intersectionHandleGroup = this.objects["intersectionHandleGroup"];
	intersectionHandleGroup.scale.set(C, C, C);
	intersectionHandleGroup.updateMatrix();

	renderHandleGroup = this.objects["renderHandleGroup"];
	renderHandleGroup.scale.set(C, C, C);
	renderHandleGroup.updateMatrix();

	group.updateMatrix();

	this.glRenderer.render(this.glScene, this.glCamera);
}

// get position of emission point
LaserPointer.prototype.getPoint = function()
{
	var pLocal = new THREE.Vector3();
	pLocal.copy(this.objects['proxy'].position);
	var pWorld = pLocal.applyMatrix4( this.group.matrix );
	return pWorld;
}


// set world position of group origin
LaserPointer.prototype.setPosition = function(position)
{
	var group = this.objects["group"];
	group.position.copy(position);
	group.updateMatrix();
}

// set world direction of laser emitted from emission point
LaserPointer.prototype.setDirection = function(direction)
{
	// rotate local y-direction of group into specified direction
	var group = this.objects["group"];
	var currentDir = this.getDirection();
	var rot = new THREE.Quaternion();
	rot.setFromUnitVectors(currentDir, direction);
	group.quaternion.multiplyQuaternions(rot, group.quaternion);
	group.updateMatrix();
}


// get direction of laser emitted from emission point
LaserPointer.prototype.getDirection = function()
{
	return this.getY();
}

LaserPointer.prototype.getX = function()
{
	var worldX = new THREE.Vector3();
	worldX.copy(this.X);
	worldX.applyQuaternion(this.group.quaternion);
	return worldX;
}

LaserPointer.prototype.getY = function()
{
	var worldY = new THREE.Vector3();
	worldY.copy(this.Y);
	worldY.applyQuaternion(this.group.quaternion);
	return worldY;
}

LaserPointer.prototype.getZ = function()
{
	var worldZ = new THREE.Vector3();
	worldZ.copy(this.Z);
	worldZ.applyQuaternion(this.group.quaternion);
	return worldZ;
}


LaserPointer.prototype.onMouseMove = function(event)
{
	// calculate mouse position in normalized device coordinates
	// (-1 to +1) for both components
	this.mouse.x =   (event.clientX / window.innerWidth)*2 - 1;
	this.mouse.y = - (event.clientY / window.innerHeight)*2 + 1;

	var mouseShift = new THREE.Vector2();
	mouseShift.copy(this.mouse).sub(this.mousePrev);
	console.log('\nmouseShift');
	console.log(event.clientX);
	console.log(event.clientY);
	console.log(this.mouse.x);
	console.log(this.mouse.y);
	console.log(this.mousePrev.x);
	console.log(this.mousePrev.y);
	console.log(mouseShift.x);
	console.log(mouseShift.y);

	this.mousePrev.copy(this.mouse);

	// update the picking ray with the camera and mouse position
	this.raycaster.setFromCamera( this.mouse, this.camera );

	obj = this.objects;
	group = obj["group"];

	var camDist = new THREE.Vector3();

	if ( this.SELECTED )
	{
		var backplaneIntersection = this.raycaster.intersectObject( obj['backPlane'] );
		if ( backplaneIntersection.length > 0 )
		{
			var planeHitpoint = backplaneIntersection[0].point;

			var shiftRelative = new THREE.Vector3(); // relative to mouse-down hit
			shiftRelative.copy(planeHitpoint).sub(this.mousedownOffset);

			var shiftAbsolute = new THREE.Vector3(); // resulting group translation in plane 
			shiftAbsolute.copy(shiftRelative).sub(group.position);
		
			if (this.SELECTED == obj['xHandleIntersection'])
			{
				var moveX = shiftAbsolute.dot(this.getX());
				var xTranslation = new THREE.Vector3();
				xTranslation.copy(this.getX()).multiplyScalar(moveX);
				group.position.add(xTranslation);
			}
			else if (this.SELECTED == obj['yHandleIntersection'])
			{
				var moveY = shiftAbsolute.dot(this.getY());
				var yTranslation = new THREE.Vector3();
				yTranslation.copy(this.getY()).multiplyScalar(moveY);
				group.position.add(yTranslation);
			}
			else if (this.SELECTED == obj['zHandleIntersection'])
			{
				var moveZ = shiftAbsolute.dot(this.getZ());
				var zTranslation = new THREE.Vector3();
				zTranslation.copy(this.getZ()).multiplyScalar(moveZ);
				group.position.add(zTranslation);
			}

			else if (this.SELECTED == obj['xRotHandleIntersection'])
			{
				var rotAngle = 0.5 * Math.PI * mouseShift.length();
				var rotX = new THREE.Quaternion();
				rotX.setFromAxisAngle(this.getX(), rotAngle);
				group.quaternion.multiplyQuaternions(rotX, group.quaternion);
			}
			else if (this.SELECTED == obj['zRotHandleIntersection'])
			{
				console.log('\nzRot');
				console.log(mouseShift.x);
				console.log(mouseShift.y);

				var rotAngle = 0.5 * Math.PI * mouseShift.length();
				console.log(rotAngle);

				var rotZ = new THREE.Quaternion();
				rotZ.setFromAxisAngle(this.getZ(), rotAngle);
				group.quaternion.multiplyQuaternions(rotZ, group.quaternion);
			}

			group.updateMatrix();
		}

		return true;
	}

	var intersections = this.raycaster.intersectObjects(this.intersectionHandleObjects);
	if ( intersections.length > 0 )
	{
		if ( this.INTERSECTED != intersections[0].object ) 
		{
			// We moused over a new handle: update backplane to prepare for possible selection:
			group = obj['group'];
			obj['backPlane'].lookAt( this.camera.position );
			obj['backPlane'].position.copy( group.position );
		}
	}

	return false;
}

// TODO:

//   beginnings of pathtracer
//      - initially just one bounce per pixel rays, to visualize distance field

//   fix laser pointer issues

//   set up webGL dropdown UI thing, use initially e.g. for scene selection

//   


LaserPointer.prototype.onMouseDown = function(event)
{
	this.mousePrev.x =   (event.clientX / window.innerWidth)*2 - 1;
	this.mousePrev.y = - (event.clientY / window.innerHeight)*2 + 1;

	// update the picking ray with the camera and mouse position
	this.raycaster.setFromCamera( this.mouse, this.camera );

	obj = this.objects;
	var intersections = this.raycaster.intersectObjects(this.intersectionHandleObjects);

	if ( intersections.length > 0 )
	{
		controls.enabled = false;
		this.SELECTED = intersections[0].object;

		// Indicate selection by making corresponding render handle slightly emissive
		var rname = this.intersectionHandleNameToRenderHandleName[this.SELECTED.name];
		var robj = this.objects[rname];
		robj.material.emissive.set( 0x404040 );

		for (var n=0; n<this.renderHandleObjects.length; n++) 
		{
			var robj_prime = this.renderHandleObjects[n];
			if (robj_prime != robj) 
			{
				robj_prime.material.emissive.set( 0x000000 );
			}
		}

		var backplaneIntersection = this.raycaster.intersectObject(obj['backPlane']);
		if ( backplaneIntersection.length > 0 ) 
		{
			var planeHitpoint = backplaneIntersection[0].point;

			// Record relative offset of mouse intersection and camera-aligned backplane
			// (which we can assume is fixed in screen space during an active manipulation)
			this.mousedownOffset.copy(planeHitpoint).sub(obj['backPlane'].position);

			group = obj['group'];
			this.groupPositionOnMouseDown   = group.position;
			this.groupQuaternionOnMouseDown = group.quaternion;
		}
	}
}

LaserPointer.prototype.onMouseUp = function(event)
{
	controls.enabled = true;
	this.SELECTED = null;


}



//////////////////////////////////////////////////////////////////////
// RayState
//////////////////////////////////////////////////////////////////////

var RayState = function(size) 
{
	this.size = size;

	var posData = new Float32Array(size*size*4); // ray position
	var dirData = new Float32Array(size*size*4); // ray direction
	var rngData = new Float32Array(size*size*4); // Random number seed
	var rgbData = new Float32Array(size*size*4); // Ray color, and wavelength

	for (var i = 0; i<size*size; ++i)
	{
		dirData[i*4 + 0] = 1.0;
		dirData[i*4 + 1] = 0.0;
		dirData[i*4 + 2] = 0.0;
		dirData[i*4 + 3] = 0.0;
		for (var t = 0; t<4; ++t)
		{
			rgbData[i*4 + t] = Math.random();
			rngData[i*4 + t] = Math.random()*4194167.0;
		}
	}

	this.posTex = new GLU.Texture(size, size, 4, true, false, true, posData);
	this.dirTex = new GLU.Texture(size, size, 4, true, false, true, dirData);
	this.rngTex = new GLU.Texture(size, size, 4, true, false, true, rngData);
	this.rgbTex = new GLU.Texture(size, size, 4, true, false, true, rgbData);
}


RayState.prototype.bind = function(shader)
{
	this.posTex.bind(0);
	this.dirTex.bind(1);
	this.rngTex.bind(2);
	this.rgbTex.bind(3);
	shader.uniformTexture("PosData", this.posTex);
	shader.uniformTexture("DirData", this.dirTex);
	shader.uniformTexture("RngData", this.rngTex);
	shader.uniformTexture("RgbData", this.rgbTex);
}


RayState.prototype.attach = function(fbo)
{
	var gl = GLU.gl;
	fbo.attachTexture(this.posTex, 0);
	fbo.attachTexture(this.dirTex, 1);
	fbo.attachTexture(this.rngTex, 2);
	fbo.attachTexture(this.rgbTex, 3);
	if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) 
	{
		GLU.fail("Invalid framebuffer");
	}
}


RayState.prototype.detach = function(fbo)
{
	var gl = GLU.gl;
	fbo.detachTexture(0);
	fbo.detachTexture(1);
	fbo.detachTexture(2);
	fbo.detachTexture(3);
}


//////////////////////////////////////////////////////////////////////
// Renderer
//////////////////////////////////////////////////////////////////////

var Renderer = function()
{
	this.gl = GLU.gl;
	var gl = GLU.gl;

	var render_canvas = document.getElementById('render-canvas');
	render_canvas.width  = window.innerWidth;
	render_canvas.height = window.innerHeight;
	this.width = render_canvas.width;
	this.height = render_canvas.height;
	
	this.initialized = false;

	// Initialize THREE.js camera
	{
		var VIEW_ANGLE = 45;
		var ASPECT = this.width / this.height ;
		var NEAR = 0.05;
		var FAR = 1000;
		this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
		this.camera.position.z = 50;
		this.camera.updateProjectionMatrix();
	}

	////////////////////////////////////////////////////////////
	// Setup three.js GL renderer
	////////////////////////////////////////////////////////////

	var ui_canvas = document.getElementById('ui-canvas');
	ui_canvas.style.top = 0;
	ui_canvas.style.position = 'absolute' 

	this.glRenderer = new THREE.WebGLRenderer( { canvas: ui_canvas,
											     alpha: true,
											     antialias: true } );
	this.glRenderer.setClearColor( 0x000000, 0 ); // the default
	this.glRenderer.setSize(this.width, this.height);
	this.glScene = new THREE.Scene();
	this.glScene.add(this.camera);

	var pointLight = new THREE.PointLight(0xFFFFFF);
	pointLight.position.x = 10;
	pointLight.position.y = 50;
	pointLight.position.z = 130;
	this.glScene.add(pointLight);

	document.body.appendChild(this.glRenderer.domElement);

	////////////////////////////////////////////////////////////
	// Setup Laser pointer
	////////////////////////////////////////////////////////////
	this.laser = new LaserPointer(this.glRenderer, this.glScene, this.camera);
	this.laser.setPosition(new THREE.Vector3(-10.0, 0.0, 0.0));
	this.laser.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));

	// @todo: define 'scenes', which consist of:

	//   - a chunk of glsl defining an SDF
	//   - starting point for laser pointer and camera

	// User settings:
	//   - switch scene
	//   - change laser pointer spread
	//   - change spectrum of emission
	//   - change global IOR?

	// @todo: some 'debug display' of the object geometry, via distance tracing.


	// @todo  for version 2.0
	//  -- more general scenes with:
	//       - a different microfacet BSDF on each object surface (both reflection and transmission)
	//       - different (but constant) IOR in SDF<0 interiors, per 'object'
	//	-- ability to switch to full pathtracing mode (bidirectional)
	


	////////////////////////////////////////////////////////////
	// Setup keypress events
	////////////////////////////////////////////////////////////
	document.body.addEventListener("keydown", this, false);
	this.glRenderer.domElement.addEventListener( 'mousemove', this, false );
	this.glRenderer.domElement.addEventListener( 'mousedown', this, false );
	this.glRenderer.domElement.addEventListener( 'mouseup',   this, false );

	this.container = document.getElementById('container');
	{
		this.stats = new Stats();
		this.stats.domElement.style.position = 'absolute';
		this.stats.domElement.style.top = '0px';
		this.container.appendChild( this.stats.domElement );
	}

	///////////////////////////////////////
	// Read shaders
	///////////////////////////////////////

	var renderer = this;
	var shaderSources = GLU.resolveShaderSource(["init", "trace", "line", "comp", "pass"]);
	if (!renderer.initialized)
	{
		renderer.setup(shaderSources);
		renderer.initialized = true;
	}
}

Renderer.prototype.handleEvent = function(event)
{
	switch (event.type)
	{
		case 'mousemove': this.onDocumentMouseMove(event); break;
		case 'mousedown': this.onDocumentMouseDown(event); break;
		case 'mouseup':   this.onDocumentMouseUp(event);   break;
		case 'keydown':   this.onKeydown(event);           break;
	}
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


Renderer.prototype.reset = function()
{
	if (!this.needsReset)
	    return;

	this.needsReset = false;
	this.wavesTraced = 0;
	this.raysTraced = 0;
	this.samplesTraced = 0;
	this.pathLength = 0;

	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.screenBuffer, 0);
	this.gl.clear(this.gl.COLOR_BUFFER_BIT);
	this.fbo.unbind();
}


Renderer.prototype.resetActiveBlock = function()
{
	//this.activeBlock = 4;
	this.activeBlock = this.raySize;
}

Renderer.prototype.setup = function(shaderSources)
{
	var gl = GLU.gl;

	// Quad VBO for rendering textures
	this.quadVbo = this.createQuadVbo();

	// shaderSources is a dict from name (e.g. "viewport")
	// to a dict {v:vertexShaderSource, f:fragmentShaderSource}
	this.initProgram  = new GLU.Shader('init',  shaderSources);
	this.traceProgram = new GLU.Shader('trace', shaderSources);
	this.lineProgram  = new GLU.Shader('line',  shaderSources);
	this.compProgram  = new GLU.Shader('comp',  shaderSources);
	this.passProgram  = new GLU.Shader('pass',  shaderSources);
	
	// table of 256 vec4 RGB colors, corresponding to the 256 wavelength samples
	// between 360.0 and 750.0 nanometres
	// @todo:  for now we will just assume a flat emission spectrum, i.e pure white light.
	this.LAMBDA_MIN = 360.0;
    this.LAMBDA_MAX = 750.0;
	this.spectrumTable = wavelengthToRgbTable();
	this.spectrum = new GLU.Texture(this.spectrumTable.length/4, 1, 4, true,  true, true, this.spectrumTable);
  
	//gl.viewport(0, 0, this.width, this.height);
	this.raySize = 32;
	this.resetActiveBlock();
	this.rayCount = this.raySize*this.raySize;
	this.currentState = 0;
	this.needsReset = true;
	this.maxPathLength = 32;
	this.rayStates = [new RayState(this.raySize), new RayState(this.raySize)];
		
	// Create the buffer of texture coordinates, which maps each drawn line
	// to its corresponding texture lookup.
	{
		this.rayVbo = new GLU.VertexBuffer();
		this.rayVbo.addAttribute("TexCoord", 3, gl.FLOAT, false);
		this.rayVbo.init(this.rayCount*2);
		var vboData = new Float32Array(this.rayCount*2*3);
		for (var i=0; i<this.rayCount; ++i)
		{
			var u = ((i % this.raySize) + 0.5) / this.raySize;
			var v = (Math.floor(i/this.raySize) + 0.5) / this.raySize;
			vboData[i*6 + 0] = vboData[i*6 + 3] = u;
			vboData[i*6 + 1] = vboData[i*6 + 4] = v;
			vboData[i*6 + 2] = 0.0;
			vboData[i*6 + 5] = 1.0;
		}
		this.rayVbo.copy(vboData);
	}

	this.fbo = new GLU.RenderTarget();

	gl.clearColor(0.0, 0.0, 0.0, 1.0);
	gl.blendFunc(gl.ONE, gl.ONE);

	this.resize(this.width, this.height);
}


Renderer.prototype.getCamera = function()
{
	return this.camera;
}


Renderer.prototype.resize = function(width, height) 
{
	this.width = width;
	this.height = height;
	this.camera.aspect = width / height;
	this.camera.updateProjectionMatrix();
	this.screenBuffer = new GLU.Texture(this.width, this.height, 4, true, false, true, null);
	this.waveBuffer   = new GLU.Texture(this.width, this.height, 4, true, false, true, null);
	this.resetActiveBlock();
	this.reset();

	this.glRenderer.setSize(width, height);
}


Renderer.prototype.composite = function()
{
	this.screenBuffer.bind(0);
	this.compProgram.bind();
	this.quadVbo.bind();
	this.compProgram.uniformTexture("Frame", this.screenBuffer);

	// Tonemap to effectively divide by total number of emitted photons
	// (and also apply gamma correction)
	this.compProgram.uniformF("Exposure", 100.0/(Math.max(this.samplesTraced, 1)));//this.raySize*this.activeBlock)));
	
	this.quadVbo.draw(this.compProgram, this.gl.TRIANGLE_FAN);
}


// With y-axis as the polar north
Renderer.prototype.sphericalPolar = function(theta, phi)
{
	Sp = Math.sin(phi);
	Cp = Math.cos(phi);
	St = Math.sin(theta);
	Ct = Math.cos(theta);
	return new THREE.Vector3( Cp*St, Ct, Sp*St );
}


Renderer.prototype.render = function()
{
	if (!this.initialized) return;

	var gl = this.gl;

	this.needsReset = true;

	var current = this.currentState;
	var next    = 1 - current;

	this.fbo.bind();

	// trace rays
	gl.viewport(0, 0, this.raySize, this.raySize);
	//gl.scissor(0, 0, this.raySize, this.activeBlock);
	//gl.enable(gl.SCISSOR_TEST);

	// We will write the next state's ray data
	this.fbo.drawBuffers(4);

	this.rayStates[next].attach(this.fbo);

	this.quadVbo.bind();

	// initialize emitter rays (the beginning of a 'wave')
	if (this.pathLength == 0)
	{
		// Start all rays at emission point(s)
		this.initProgram.bind();

		// Read random seed from the current state
		this.rayStates[current].rngTex.bind(0); 
		this.initProgram.uniformTexture("RngData", this.rayStates[current].rngTex);

		// Read wavelength -> RGB table 
		this.spectrum.bind(1);
		this.initProgram.uniformTexture("Spectrum", this.spectrum);
		
		// Emitter data (currently just location and direction)
		emitterPos = this.laser.getPoint();
		this.initProgram.uniform3F("EmitterPos", emitterPos.x, emitterPos.y, emitterPos.z);

		emitterDir = this.laser.getDirection();
		this.initProgram.uniform3F("EmitterDir", emitterDir.x, emitterDir.y, emitterDir.z);

		// Write emitted ray initial conditions into 'next' state
		this.quadVbo.draw(this.initProgram, gl.TRIANGLE_FAN);

		// Make this initial state be 'current'
		current = 1 - current;
		next    = 1 - next;

		// And we prepare to write into the 'next' state
		this.fbo.drawBuffers(4);

		this.rayStates[next].attach(this.fbo);
	}

	// Raytrace into scene, generating new ray pos/dir data in 'next' rayStates textures
	{
		this.traceProgram.bind();

		// Use the current state as the initial conditions
		this.rayStates[current].bind(this.traceProgram);

		// Generate the next ray state
		this.quadVbo.draw(this.traceProgram, gl.TRIANGLE_FAN);

		this.rayStates[next].detach(this.fbo);
	}

	// Draw the next set of lines into the wave buffer
	{
		//gl.disable(gl.SCISSOR_TEST);
		gl.viewport(0, 0, this.width, this.height);

		this.fbo.drawBuffers(1);
		this.fbo.attachTexture(this.waveBuffer, 0);

		if (this.pathLength == 0 || this.wavesTraced==0)
		{
			// Clear wavebuffer before the first bounce
			gl.clear(gl.COLOR_BUFFER_BIT);
		}

		gl.enable(gl.BLEND);

		this.lineProgram.bind();

		// Setup projection matrix
		var projectionMatrix = this.camera.projectionMatrix.toArray();
		var projectionMatrixLocation = this.lineProgram.getUniformLocation("u_projectionMatrix");
		gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

		// Setup modelview matrix (to match camera)
		this.camera.updateMatrixWorld();
		var matrixWorldInverse = new THREE.Matrix4();
		matrixWorldInverse.getInverse( this.camera.matrixWorld );
		var modelViewMatrix = matrixWorldInverse.toArray();
		var modelViewMatrixLocation = this.lineProgram.getUniformLocation("u_modelViewMatrix");
		gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);


		this.rayStates[current].posTex.bind(0); // PosDataA = current.posTex
		this.rayStates[   next].posTex.bind(1); // PosDataB = next.posTex
		this.rayStates[current].rgbTex.bind(2); // current  = current.rgbTex

		this.lineProgram.uniformTexture("PosDataA", this.rayStates[current].posTex);
		this.lineProgram.uniformTexture("PosDataB", this.rayStates[   next].posTex);
		this.lineProgram.uniformTexture("RgbData",  this.rayStates[current].rgbTex);

		this.rayVbo.bind(); // Binds the TexCoord attribute
		this.rayVbo.draw(this.lineProgram, gl.LINES, this.raySize*this.activeBlock*2);

		this.raysTraced += this.raySize*this.activeBlock;
		this.pathLength += 1;
	}

	this.quadVbo.bind();

	// Update the screenBuffer with the waveBuffer contents
	if (this.pathLength==this.maxPathLength || this.wavesTraced==0)
	{
		this.fbo.attachTexture(this.screenBuffer, 0);
		this.waveBuffer.bind(0);
		this.passProgram.bind();
		this.passProgram.uniformTexture("Frame", this.waveBuffer);
		this.quadVbo.draw(this.passProgram, gl.TRIANGLE_FAN);

		this.samplesTraced += this.raySize*this.activeBlock;

		if (this.pathLength == this.maxPathLength)
		{
			this.wavesTraced += 1;
			this.pathLength = 0;
		}
	}

	gl.disable(gl.BLEND);
	
	this.fbo.unbind();

	// Final composite of screenBuffer to window
	this.composite();

	// Render laser pointer
	this.laser.render(this.glScene);

	// Update raytracing state
	this.currentState = next;

	this.stats.update();
}


Renderer.prototype.onDocumentMouseMove = function(event)
{
	if (!this.initialized) return;
	event.preventDefault();
	if (this.laser.onMouseMove(event)) this.reset();
}

Renderer.prototype.onDocumentMouseDown = function(event)
{
	if (!this.initialized) return;
	event.preventDefault();
	this.laser.onMouseDown(event);
}

Renderer.prototype.onDocumentMouseUp = function(event)
{
	if (!this.initialized) return;
	event.preventDefault();
	this.laser.onMouseUp(event);
}

Renderer.prototype.onKeydown = function(event)
{
	if (!this.initialized) return;
	event.preventDefault();

	var charCode = (event.which) ? event.which : event.keyCode;
	var fKeyCode = 70;

	if (charCode == fKeyCode)
	{
		var element	= document.body;
		if ( 'webkitCancelFullScreen' in document )
		{
			element.webkitRequestFullScreen();
		}
		else if ( 'mozCancelFullScreen' in document )
		{
			element.mozRequestFullScreen();
		}
		else
		{
			console.assert(false);
		}
	}
}



