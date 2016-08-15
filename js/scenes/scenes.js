

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
	var sceneObj = snelly.getLoadedScene();
	var sceneScale = sceneObj.getScale();
	laser.setEmissionRadius(0.01*sceneScale);
	laser.setEmissionSpreadAngle(5.0);
	laser.buildEmitterGeo();
}






