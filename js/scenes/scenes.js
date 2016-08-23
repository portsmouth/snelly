

function Scene(name, desc)
{
	this._name = name;
	this._desc = desc;
	this._settings = {};
}

Scene.prototype.getName = function()
{
	return this._name;
}

Scene.prototype.getDesc = function()
{
	return this._desc;
}

Scene.prototype.getSettings = function()
{
	return this._settings;
}

// Initial emitter defaults
Scene.prototype.setLaser = function(laser)
{
	laser.buildEmitterGeo();
}

Scene.prototype.getBox = function()
{
	var s = 2.0*this.getScale();
	var min = new THREE.Vector3(-s, -s, -s);
	var max = new THREE.Vector3( s,  s,  s);
	return new THREE.Box3(min, max);
}






