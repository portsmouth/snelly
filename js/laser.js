

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
	this.mousedownOffset = new THREE.Vector3();

	this.INTERSECTED = null;
	this.SELECTED = null;

	this.objects = {};

	// Initialize housing geometry:
	this.setEmissionRadius(0.1);
	this.setEmissionSpreadAngle(5.0);
	this.buildHousing();

	///////////////////////////////////////////////////////////
	// Intersection geometry for handle controls (invisible)
	///////////////////////////////////////////////////////////

	this.intersectionHandleGroup = new THREE.Object3D();
	var intersectionHandleGroup = this.intersectionHandleGroup;

	var rI = 0.4;
	var hI = 8.0;
	var RI = 2.0;

	// x-handle
	var xHandleIntersectionGeo       = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var xHandleIntersectionMaterial  = new THREE.MeshPhongMaterial( { color: 'red', visible: false } );
	var xHandleIntersectionObj       = new THREE.Mesh( xHandleIntersectionGeo, xHandleIntersectionMaterial  );
	xHandleIntersectionObj.rotation.z = 0.5*Math.PI;
	xHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(xHandleIntersectionObj);

	// y-handle
	var yHandleIntersectionGeo        = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var yHandleIntersectionMaterial   = new THREE.MeshPhongMaterial( { color: 'green', visible: false } );
	var yHandleIntersectionObj        = new THREE.Mesh( yHandleIntersectionGeo, yHandleIntersectionMaterial  );
	yHandleIntersectionObj.rotation.y = 0.5*Math.PI;
	yHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(yHandleIntersectionObj);

	// z-handle
	var zHandleIntersectionGeo        = new THREE.CylinderGeometry( rI, rI, hI, 32 );
	var zHandleIntersectionMaterial   = new THREE.MeshPhongMaterial( { color: 'blue', visible: false } );
	var zHandleIntersectionObj        = new THREE.Mesh( zHandleIntersectionGeo, zHandleIntersectionMaterial  );
	zHandleIntersectionObj.rotation.x = 0.5*Math.PI;
	zHandleIntersectionObj.position.x = -0.0;
	intersectionHandleGroup.add(zHandleIntersectionObj);

	// x-rot-handle
	var xRotHandleIntersectionGeo      = new THREE.TorusGeometry( RI, rI, 32, 100 );
	var xRotHandleIntersectionMaterial = new THREE.MeshPhongMaterial( { color: 'red', visible: false } );
	var xRotHandleIntersectionObj      = new THREE.Mesh( xRotHandleIntersectionGeo, xRotHandleIntersectionMaterial  );
	xRotHandleIntersectionObj.rotation.y = 0.5*Math.PI;
	intersectionHandleGroup.add(xRotHandleIntersectionObj);

	// z-rot-handle
	var zRotHandleIntersectionGeo      = new THREE.TorusGeometry( RI, rI, 32, 100 );
	var zRotHandleIntersectionMaterial = new THREE.MeshPhongMaterial( { color: 'blue', visible: false } );
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
	this.objects['yHandleIntersection'] = xHandleIntersectionObj;
	this.objects['zHandleIntersection'] = zHandleIntersectionObj;
	this.objects['xRotHandleIntersection'] = xRotHandleIntersectionObj;
	this.objects['zRotHandleIntersection'] = zRotHandleIntersectionObj;
	this.objects['intersectionHandleGroup'] = intersectionHandleGroup;

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

	// Initialize position and direction
	this.setDirection(new THREE.Vector3(1.0, 0.0, 0.0));
	this.setPosition(new THREE.Vector3(-10.0, 0.0, 0.0));
}


LaserPointer.prototype.buildHousing = function()
{
	var group = this.group;
	if ('housing' in this.objects)
	{
		var housing = this.objects['housing'];
		group.remove(housing);
	}
	
	// Laser housing geometry
	var rPointer = 3.0*this.getEmissionRadius();
	var lPointer = 10.0*rPointer;
	var housingGeo      = new THREE.CylinderGeometry( rPointer, rPointer, lPointer, 32 );
	var housingMaterial = new THREE.MeshPhongMaterial( { color: 0xdddddd, specular: 0x999999, shininess: 60 } );
	var housingObj      = new THREE.Mesh( housingGeo, housingMaterial  );
	group.add(housingObj);
	this.housingObj = housingObj;
	this.objects['housing'] = housingObj;

	// emission point proxy
	var proxyGeo      = new THREE.SphereGeometry(this.getEmissionRadius());
	var proxyMaterial = new THREE.MeshPhongMaterial( { color: 0xff0000, visible: 'true' } );
	var proxyObj      = new THREE.Mesh( proxyGeo, proxyMaterial  );
	proxyObj.position.y = 0.5*lPointer;
	group.add(proxyObj);
	this.objects['proxy'] = proxyObj;
}

LaserPointer.prototype.render = function()
{
	// Ensure laser handles are scaled according to current camera position
	// (to be roughly constant size in screen space)
	var camDist = new THREE.Vector3();
	camDist.copy(this.objects["group"].position).sub(this.camera.position);
	var C = 0.05*camDist.length();

	intersectionHandleGroup = this.objects["intersectionHandleGroup"];
	intersectionHandleGroup.scale.set(C, C, C);
	intersectionHandleGroup.updateMatrix();

	renderHandleGroup = this.objects["renderHandleGroup"];
	renderHandleGroup.scale.set(C, C, C);
	renderHandleGroup.updateMatrix();

	this.glRenderer.render(this.glScene, this.glCamera);
}

/// Queries:

// get position of emission point
LaserPointer.prototype.getPoint = function()
{
	var pLocal = new THREE.Vector3();
	pLocal.copy(this.objects['proxy'].position);
	var pWorld = pLocal.applyMatrix4( this.group.matrix );
	return pWorld;
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

LaserPointer.prototype.setEmissionRadius = function(radius)
{
	this.emissionRadius = radius;
	this.buildHousing();
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

	//this.mouse.x =   (event.clientX / window.innerWidth)*2 - 1;
	//this.mouse.y = - (event.clientY / window.innerHeight)*2 + 1;

	var mouseShift = new THREE.Vector2();
	mouseShift.copy(this.mouse).sub(this.mousePrev);
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

			// @todo: make ring always rotate in the direction
			//        it is being 'dragged' (could just flip depending
			//        on camDir.ringNormal)

			else if (this.SELECTED == obj['xRotHandleIntersection'])
			{
				var rotAngle = 0.5 * Math.PI * (-mouseShift.x -mouseShift.y);
				var rotX = new THREE.Quaternion();
				rotX.setFromAxisAngle(this.getX(), rotAngle);
				group.quaternion.multiplyQuaternions(rotX, group.quaternion);
			}
			else if (this.SELECTED == obj['zRotHandleIntersection'])
			{
				var rotAngle = 0.5 * Math.PI * (-mouseShift.x -mouseShift.y);
				var rotZ = new THREE.Quaternion();
				rotZ.setFromAxisAngle(this.getZ(), rotAngle);
				group.quaternion.multiplyQuaternions(rotZ, group.quaternion);
			}

			group.updateMatrix();
		}

		return true;
	}

	for (var n=0; n<this.renderHandleObjects.length; n++) 
	{
		this.renderHandleObjects[n].material.emissive.set( 0x000000 );
		this.intersectionHandleObjects[n].material.emissive.set( 0x000000 );
	}

	var intersections = this.raycaster.intersectObjects(this.intersectionHandleObjects);
	if ( intersections.length > 0 )
	{
		var intersected = intersections[0].object;
		var rname = this.intersectionHandleNameToRenderHandleName[intersected.name];

		this.objects[intersected.name].material.emissive.set( 0x404040 );
		this.objects[rname].material.emissive.set( 0x404040 );

		if ( this.INTERSECTED != intersected ) 
		{
			// We moused over a new handle: update backplane to prepare for possible selection:
			group = obj['group'];
			obj['backPlane'].lookAt( this.camera.position );
			obj['backPlane'].position.copy( group.position );
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
	this.controls.enabled = true;
	this.SELECTED = null;
}


