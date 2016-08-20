

//////////////////////////////////////////////////////////////////////
// Laser pointer UI
//////////////////////////////////////////////////////////////////////

var LaserPointer = function(glRenderer, glScene, glCamera, controls) 
{
	// Main group of all objects moving rigidly with the laser pointer:
	this.group = new THREE.Object3D();
	var group = this.group;

	this.camera = glCamera;
	this.raycaster = new THREE.Raycaster();
	this.controls = controls;

	this.mouse     = new THREE.Vector2();
	this.mousePrev = new THREE.Vector2();
	this.mousedownPlaneHitpoint = new THREE.Vector3();

	this.INTERSECTED = null;
	this.SELECTED = null;

	this.objects = {};

	///////////////////////////////////////////////////////////
	// Intersection geometry for handle controls (invisible)
	///////////////////////////////////////////////////////////

	this.intersectionHandleGroup = new THREE.Object3D();
	var intersectionHandleGroup = this.intersectionHandleGroup;

	var rI = 0.4;
	var hI = 12.0;
	var RI = 2.0;

	var intersectionDebug = false;

	// x-handle
	var xHandleIntersectionGeo       = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var xHandleIntersectionMaterial  = new THREE.MeshPhongMaterial( { color: 'red', visible: intersectionDebug } );
	var xHandleIntersectionObj       = new THREE.Mesh( xHandleIntersectionGeo, xHandleIntersectionMaterial  );
	xHandleIntersectionObj.rotation.z = 0.5*Math.PI;
	xHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(xHandleIntersectionObj);

	// y-handle
	var yHandleIntersectionGeo        = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var yHandleIntersectionMaterial   = new THREE.MeshPhongMaterial( { color: 'green', visible: intersectionDebug } );
	var yHandleIntersectionObj        = new THREE.Mesh( yHandleIntersectionGeo, yHandleIntersectionMaterial  );
	yHandleIntersectionObj.rotation.y = 0.5*Math.PI;
	yHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(yHandleIntersectionObj);

	// z-handle
	var zHandleIntersectionGeo        = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var zHandleIntersectionMaterial   = new THREE.MeshPhongMaterial( { color: 'blue', visible: intersectionDebug } );
	var zHandleIntersectionObj        = new THREE.Mesh( zHandleIntersectionGeo, zHandleIntersectionMaterial  );
	zHandleIntersectionObj.rotation.x = 0.5*Math.PI;
	zHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(zHandleIntersectionObj);

	// x-rot-handle
	var xRotHandleIntersectionGeo      = new THREE.TorusGeometry( RI, rI, 32, 100 );
	var xRotHandleIntersectionMaterial = new THREE.MeshPhongMaterial( { color: 'red', visible: intersectionDebug } );
	var xRotHandleIntersectionObj      = new THREE.Mesh( xRotHandleIntersectionGeo, xRotHandleIntersectionMaterial  );
	xRotHandleIntersectionObj.rotation.y = 0.5*Math.PI;
	intersectionHandleGroup.add(xRotHandleIntersectionObj);

	// z-rot-handle
	var zRotHandleIntersectionGeo      = new THREE.TorusGeometry( RI, rI, 32, 100 );
	var zRotHandleIntersectionMaterial = new THREE.MeshPhongMaterial( { color: 'blue', visible: intersectionDebug } );
	var zRotHandleIntersectionObj      = new THREE.Mesh( zRotHandleIntersectionGeo, zRotHandleIntersectionMaterial  );
	zRotHandleIntersectionObj.rotation.y = 0.0;
	intersectionHandleGroup.add(zRotHandleIntersectionObj);

	group.add(intersectionHandleGroup);

	// translater intersection geo, for dragging
	var translaterGeo      = new THREE.SphereGeometry(0.98*RI, 32, 32);
	var translaterMaterial = new THREE.MeshPhongMaterial( { color: 0xff0000, visible: true } );
	var translaterObj      = new THREE.Mesh( translaterGeo, translaterMaterial  );
	group.add(translaterObj);


	//////////////////////////////////////////////////////////////////////////////
	// Rendered geometry for handle controls (thin, visually appealing lines)
	//////////////////////////////////////////////////////////////////////////////
	
	this.renderHandleGroup = new THREE.Object3D();
	var renderHandleGroup = this.renderHandleGroup;

	var rR = 0.06;
	var hR = 8.0;
	var RR = 2.0;

	// x-handle
	var xHandleRenderGeo      = new THREE.CylinderGeometry( rR, rR, hR, 32 );
	var xHandleRenderMaterial = new THREE.MeshPhongMaterial( { color: 0x800000, shininess: 10 } );
	var xHandleRenderObj      = new THREE.Mesh( xHandleRenderGeo, xHandleRenderMaterial  );
	xHandleRenderObj.rotation.z = 0.5*Math.PI;
	xHandleRenderObj.position.x = -0.0;
	renderHandleGroup.add(xHandleRenderObj);

	// y-handle
	var yHandleRenderGeo      = new THREE.CylinderGeometry( rR, rR, hR, 32 );
	var yHandleRenderMaterial = new THREE.MeshPhongMaterial( { color: 0x008000, shininess: 10 } );
	var yHandleRenderObj      = new THREE.Mesh( yHandleRenderGeo, yHandleRenderMaterial  );
	yHandleRenderObj.rotation.y = 0.5*Math.PI;
	yHandleRenderObj.position.x = -0.0;
	renderHandleGroup.add(yHandleRenderObj);

	// z-handle
	var zHandleRenderGeo      = new THREE.CylinderGeometry( rR, rR, hR, 32 );
	var zHandleRenderMaterial = new THREE.MeshPhongMaterial( { color: 0x000080, shininess: 30 } );
	var zHandleRenderObj      = new THREE.Mesh( zHandleRenderGeo, zHandleRenderMaterial  );
	zHandleRenderObj.rotation.x = 0.5*Math.PI;
	zHandleRenderObj.position.x = -0.0;
	renderHandleGroup.add(zHandleRenderObj);

	group.add(renderHandleGroup);

	// x-rot-handle
	var xRotHandleRenderGeo      = new THREE.TorusGeometry( RR, rR, 32, 100 );
	var xRotHandleMaterial = new THREE.MeshPhongMaterial( { color: 0x800000, shininess: 10 } );
	var xRotHandleRenderObj = new THREE.Mesh( xRotHandleRenderGeo, xRotHandleMaterial  );
	xRotHandleRenderObj.rotation.y = 0.5*Math.PI;
	renderHandleGroup.add(xRotHandleRenderObj);

	// z-rot-handle
	var zRotHandleRenderGeo      = new THREE.TorusGeometry( RR, rR, 32, 100 );
	var zRotHandleMaterial = new THREE.MeshPhongMaterial( { color: 0x000080, shininess: 10 } );
	var zRotHandleRenderObj = new THREE.Mesh( zRotHandleRenderGeo, zRotHandleMaterial  );
	zRotHandleRenderObj.rotation.y = 0.0;
	renderHandleGroup.add(zRotHandleRenderObj);

	
	// back plane to implement dragging 
	var backPlaneObj =  new THREE.Mesh(
						new THREE.PlaneGeometry(2000, 2000),
						new THREE.MeshBasicMaterial( { visible: false } ) );
	backPlaneObj.position.x = 0.0;
	backPlaneObj.position.y = 0.0;

	this.X = new THREE.Vector3(1, 0, 0);
	this.Y = new THREE.Vector3(0, 1, 0);
	this.Z = new THREE.Vector3(0, 0, 1);

	this.objects['group'] = this.group;

	this.objects['xHandleIntersection'] = xHandleIntersectionObj;
	this.objects['yHandleIntersection'] = yHandleIntersectionObj;
	this.objects['zHandleIntersection'] = zHandleIntersectionObj;
	this.objects['xRotHandleIntersection'] = xRotHandleIntersectionObj;
	this.objects['zRotHandleIntersection'] = zRotHandleIntersectionObj;
	this.objects['intersectionHandleGroup'] = intersectionHandleGroup;
	this.objects['translater'] = translaterObj;

	this.objects['xHandleRender'] = xHandleRenderObj;
	this.objects['yHandleRender'] = yHandleRenderObj;
	this.objects['zHandleRender'] = zHandleRenderObj;
	this.objects['xRotHandleRender'] = xRotHandleRenderObj;
	this.objects['zRotHandleRender'] = zRotHandleRenderObj;
	this.objects['renderHandleGroup'] = renderHandleGroup;

	this.objects['backPlane'] = backPlaneObj;

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
									"zRotHandleIntersection",
									"translater" ];
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

	this.resize(window.innerWidth, window.innerHeight);
}


LaserPointer.prototype.buildEmitterGeo = function()
{
	var group = this.group;
	if (typeof this.emitterObj !== 'undefined') 
	{
		group.remove(this.emitterObj);
	}

	// Laser emission geometry
	var rPointer = 1.1*this.getEmissionRadius();
	this.emitterLength = 0.01*rPointer;

	var emitterGeo      = new THREE.CylinderGeometry(rPointer, rPointer, this.emitterLength, 32, 32);
	var emitterMaterial = new THREE.MeshPhongMaterial( { color: 0xdddddd, specular: 0x999999, shininess: 60, visible: true } );
	var emitterObj      = new THREE.Mesh( emitterGeo, emitterMaterial  );
	group.add(emitterObj);
	this.emitterObj = emitterObj;
}


LaserPointer.prototype.resize = function(width, height)
{
	// Set up depth buffer rendering
	/* @todo: disabling for now
	var gl = GLU.gl;
	this.depthTarget = null;
	if (GLU.depthTexExt != null)
	{
		this.depthTarget = new THREE.WebGLRenderTarget(width, height);
		this.depthTarget.texture.format = THREE.RGBFormat;
		this.depthTarget.texture.minFilter = THREE.NearestFilter;
		this.depthTarget.texture.magFilter = THREE.NearestFilter;
		this.depthTarget.texture.generateMipmaps = false;
		this.depthTarget.stencilBuffer = false;
		this.depthTarget.depthBuffer = true;

		this.depthTarget.depthTexture = new THREE.DepthTexture();
		this.depthTarget.depthTexture.type = THREE.UnsignedShortType;
	}
	*/
}


LaserPointer.prototype.renderBox = function()
{
	// For a basic UI giving some aid to getting a sense of the 
	// orientation of the laser in 3d space:

	// Draw the rough bbox of the scene,
	// i.e. origin-centered cube of size e.g. 10*sceneScale


	// Then draw 6 dotted lines projecting the laser center onto the
	// 6 planes of the box, giving some sense of 


	// Also draw a dotted line down the laser axis itself.

}

LaserPointer.prototype.render = function()
{
	// Ensure laser handles are scaled according to current camera position
	// (to be roughly constant size in screen space)
	var camDist = new THREE.Vector3();
	camDist.copy(this.objects["group"].position).sub(this.camera.position);
	var C = 0.015*camDist.length();

	// Ensure that handles are always > emitter geo size
	// (so on zooming in, body doesn't 'swallow' the manips).
	var C = Math.max(C, 0.05*this.getEmissionRadius());

	intersectionHandleGroup = this.objects["intersectionHandleGroup"];
	intersectionHandleGroup.scale.set(C, C, C);
	intersectionHandleGroup.updateMatrix();

	translater = this.objects['translater'];
	translater.scale.set(C, C, C);
	translater.updateMatrix();

	renderHandleGroup = this.objects["renderHandleGroup"];
	renderHandleGroup.scale.set(C, C, C);
	renderHandleGroup.updateMatrix();

	// render scene to depth texture for later depth tests
	/* @todo: disabling for now
	if (GLU.depthTexExt != null)
	{
		this.glRenderer.render(this.glScene, this.glCamera, this.depthTarget);
	}
	*/

	this.glRenderer.render(this.glScene, this.glCamera);
}

/// Queries:

LaserPointer.prototype.getDepthTarget = function()
{
	return this.depthTarget;
}

// get world position of emission point
LaserPointer.prototype.getPoint = function()
{
	var pLocal = new THREE.Vector3();
	pLocal.copy(this.emitterObj.position);
	pLocal.y += 0.5*this.emitterLength;
	var pWorld = pLocal.applyMatrix4( this.group.matrix );
	return pWorld;
}

// get world direction of laser emitted from emission point
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

LaserPointer.prototype.getEmissionRadius = function()
{
	return this.emissionRadius;
}

// In degrees
LaserPointer.prototype.getEmissionSpreadAngle = function()
{
	return this.emissionSpread;
}

/// Interactions:

LaserPointer.prototype.toggleVisibility = function(visible)
{
	var group = this.objects["group"];
	group.visible = visible;
}

LaserPointer.prototype.isVisible = function(visible)
{
	var group = this.objects["group"];
	return group.visible;
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
	direction.normalize();
	rot.setFromUnitVectors(currentDir, direction);
	group.quaternion.multiplyQuaternions(rot, group.quaternion);
	group.quaternion.normalize();
	group.updateMatrix();
}

LaserPointer.prototype.setTarget = function(target)
{
	var newDir = target.clone();
	newDir.sub(this.getPoint());
	newDir.normalize();
	this.setDirection(newDir);
}

LaserPointer.prototype.setEmissionRadius = function(radius)
{
	this.emissionRadius = radius;
	this.buildEmitterGeo();
}

// In degrees
LaserPointer.prototype.setEmissionSpreadAngle = function(spreadAngleDegrees)
{
	this.emissionSpread = spreadAngleDegrees;
}

LaserPointer.prototype.onMouseMove = function(event)
{
	// calculate mouse position in normalized device coordinates
	// (-1 to +1) for both components
	this.mouse.x =   (( event.clientX - this.glRenderer.domElement.offsetLeft) / this.glRenderer.domElement.width)*2 - 1;
	this.mouse.y = - (( event.clientY - this.glRenderer.domElement.offsetTop) / this.glRenderer.domElement.height)*2 + 1;
	var mouseShift = new THREE.Vector2();
	mouseShift.copy(this.mouse).sub(this.mousePrev);
	this.mousePrev.copy(this.mouse);

	// update the picking ray with the camera and mouse position
	this.raycaster.setFromCamera( this.mouse, this.camera );
	obj = this.objects;
	group = obj["group"];
	var sceneObj = snelly.getLoadedScene();
	var sceneScale = sceneObj.getScale();

	if ( this.SELECTED != null )
	{
		var backplaneIntersection = this.raycaster.intersectObject( obj['backPlane'] );
		if ( backplaneIntersection.length > 0 )
		{
			var planeHitpoint = backplaneIntersection[0].point;
			var shift = new THREE.Vector3(); // relative to mouse-down hit
			shift.copy(planeHitpoint).sub(this.mousedownPlaneHitpoint);
			var shiftAbsolute = new THREE.Vector3();
			shiftAbsolute.copy(shift);
			var shiftDist = shift.length();
			shift.normalize();
			this.mousedownPlaneHitpoint.copy(planeHitpoint);

			var epsilonLength = 1.0e-6*sceneScale;

			// Translation on dragging 'rods'
			if (this.SELECTED == obj['xHandleIntersection'])
			{
				var moveX = shiftDist / (shift.dot(this.getX()) + epsilonLength);
				var xTranslation = new THREE.Vector3();
				xTranslation.copy(this.getX()).multiplyScalar(moveX);
				group.position.add(xTranslation);
			}
			else if (this.SELECTED == obj['yHandleIntersection'])
			{
				var moveY = shiftDist / (shift.dot(this.getY()) + epsilonLength);
				var yTranslation = new THREE.Vector3();
				yTranslation.copy(this.getY()).multiplyScalar(moveY);
				group.position.add(yTranslation);
			}
			else if (this.SELECTED == obj['zHandleIntersection'])
			{
				var moveZ = shiftDist / (shift.dot(this.getZ()) + epsilonLength);
				var zTranslation = new THREE.Vector3();
				zTranslation.copy(this.getZ()).multiplyScalar(moveZ);
				group.position.add(zTranslation);
			}

			// rotation on dragging 'rings'
			else if (this.SELECTED == obj['xRotHandleIntersection'])
			{
				// apply angular motion induced by drag
				var radiusVector = new THREE.Vector3();
				radiusVector.copy(this.SELECTION_HITPOINT);
				radiusVector.sub(group.position);
				var radius = radiusVector.length() + epsilonLength;

				var localX = new THREE.Vector3(1.0, 0.0, 0.0);
				var localY = new THREE.Vector3(0.0, 1.0, 0.0);
				var worldX = localX.transformDirection( this.camera.matrix );
				var worldY = localY.transformDirection( this.camera.matrix );
				worldX.multiplyScalar(mouseShift.x);
				worldY.multiplyScalar(mouseShift.y);
				var worldVelocity = new THREE.Vector3();
				worldVelocity.copy(worldX);
				worldVelocity.add(worldY);

				var L = new THREE.Vector3();
				L.crossVectors( radiusVector, worldVelocity );
				L.normalize();
				var rotAngleX = (shiftDist/radius) * L.dot(this.getX());
				var rotX = new THREE.Quaternion();
				rotX.setFromAxisAngle(this.getX(), rotAngleX);
				group.quaternion.multiplyQuaternions(rotX, group.quaternion);
			}
			else if (this.SELECTED == obj['zRotHandleIntersection'])
			{
				// apply angular motion induced by drag
				var radiusVector = new THREE.Vector3();
				radiusVector.copy(this.SELECTION_HITPOINT);
				radiusVector.sub(group.position);
				var radius = radiusVector.length() + epsilonLength;

				var localX = new THREE.Vector3(1.0, 0.0, 0.0);
				var localY = new THREE.Vector3(0.0, 1.0, 0.0);
				var worldX = localX.transformDirection( this.camera.matrix );
				var worldY = localY.transformDirection( this.camera.matrix );
				worldX.multiplyScalar(mouseShift.x);
				worldY.multiplyScalar(mouseShift.y);
				var worldVelocity = new THREE.Vector3();
				worldVelocity.copy(worldX);
				worldVelocity.add(worldY);

				var L = new THREE.Vector3();
				L.crossVectors( radiusVector, worldVelocity );
				L.normalize();
				var rotAngleZ = (shiftDist/radius) * L.dot(this.getZ());
				var rotZ = new THREE.Quaternion();
				rotZ.setFromAxisAngle(this.getZ(), rotAngleZ);
				group.quaternion.multiplyQuaternions(rotZ, group.quaternion);
			}

			// translation on dragging laser 'housing'
			else if (this.SELECTED == obj['translater'])
			{
				// choose 'the plane most orthogonal to the view dir'
				var camDir = this.camera.getWorldDirection();
				var xproj = Math.abs(camDir.dot(this.getX())); // votes for yz plane
				var yproj = Math.abs(camDir.dot(this.getY())); // votes for xy plane
				var zproj = Math.abs(camDir.dot(this.getZ())); // votes for xy plane
				var planeNormal;
				if (xproj>yproj)
				{
					if (xproj>zproj) planeNormal = this.getX();
					else             planeNormal = this.getZ();
				}
				else
				{
					if (yproj>zproj) planeNormal = this.getY();
					else             planeNormal = this.getZ();
				}

				// project shiftRelative on chosen plane
				var projectionDelta = new THREE.Vector3();
				projectionDelta.copy(planeNormal);

				projectionDelta.multiplyScalar(shiftAbsolute.dot(planeNormal));
				var translation = new THREE.Vector3();
				translation.copy(shiftAbsolute);
				translation.sub(projectionDelta);
				group.position.add(translation);
			}

			group.updateMatrix();
		}

		return true;
	}

	else
	{
		for (var n=0; n<this.renderHandleObjects.length; n++) 
		{
			this.renderHandleObjects[n].material.emissive.set( 0x000000 );
		}
		this.emitterObj.material.emissive.set( 0x000000 );
		translater.material.emissive.set( 0x000000 );

		var intersections = this.raycaster.intersectObjects(this.intersectionHandleObjects);
		if ( intersections.length > 0 )
		{
			var intersected = intersections[0].object;

			if (intersected == this.objects['translater'])
			{
				this.emitterObj.material.emissive.set( 0x404040 );
				translater.material.emissive.set( 0x102010 );
			}
			else
			{
				var rname = this.intersectionHandleNameToRenderHandleName[intersected.name];
				this.objects[rname].material.emissive.set( 0x404040 );
			}

			if ( this.INTERSECTED != intersected ) 
			{
				// We moused over a new handle: update backplane to prepare for possible selection:
				group = obj['group'];
				obj['backPlane'].lookAt( this.camera.position );
				obj['backPlane'].position.copy( group.position );
			}
		}
	}

	return false;
}


LaserPointer.prototype.onMouseDown = function(event)
{
	this.mousePrev.x =   (( event.clientX - this.glRenderer.domElement.offsetLeft ) / this.glRenderer.domElement.width)*2 - 1;
	this.mousePrev.y = - (( event.clientY - this.glRenderer.domElement.offsetTop ) / this.glRenderer.domElement.height)*2 + 1;

	// update the picking ray with the camera and mouse position
	this.raycaster.setFromCamera( this.mouse, this.camera );

	obj = this.objects;
	var intersections = this.raycaster.intersectObjects(this.intersectionHandleObjects);

	if ( intersections.length > 0 )
	{
		this.controls.enabled = false;
		this.SELECTED           = intersections[0].object;
		this.SELECTION_HITPOINT = intersections[0].point;
		var backplaneIntersection = this.raycaster.intersectObject(obj['backPlane']);
		if ( backplaneIntersection.length > 0 ) 
		{
			var planeHitpoint = backplaneIntersection[0].point;

			// Record relative offset of mouse intersection and camera-aligned backplane
			// (which we can assume is fixed in screen space during an active manipulation)
			this.mousedownPlaneHitpoint.copy(planeHitpoint);
			group = obj['group'];
			this.groupPositionOnMouseDown   = group.position;
			this.groupQuaternionOnMouseDown = group.quaternion;
		}
	}
}

LaserPointer.prototype.onMouseUp = function(event)
{
	this.controls.enabled = true;
	this.SELECTED = null;
}


