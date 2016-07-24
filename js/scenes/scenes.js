

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



////////////////////////////////////////////////////////////////////////
// Instantiate a variety of scenes
////////////////////////////////////////////////////////////////////////

function createScenes(renderer)
{
	renderer.addScene(new SphereScene("sphere", "Simple sphere"));
	renderer.addScene(new FibreScene("fibre", "Simple optical fibre"));



}




