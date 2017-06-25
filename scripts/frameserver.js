
// Movie generation process:

// 1) Add the following logic in pre/postframeCallback, which will 
// issue frames to the server once a sample theshold has been reached per-frame
/*
Scene.prototype.preframeCallback = function(snelly, gl)
{
	let renderer  = snelly.getRenderer();
	let camera    = snelly.getCamera();
	let controls  = snelly.getControls();
	let materials = snelly.getMaterials();
	let gui       = snelly.getGUI();

	let FPS = 24.0;
	let time = this.animFrame/FPS;
	this.endFrame = 120.0; * FPS; // For example, for a 2 minute sequence

	// animate camera (here a simple 'turntable' orbit about the original cam target)
	let axis = camera.up;

	if (this.animFrame > this.endFrame) this.animFrame = 0;
	if (this.animFrame==0)
	{
		this.advanceFrame = false;
		///
		// Do any one-time initial setup of the scene state here
		///
	}

	///
	// Animate scene state according to current time value here
	// (e.g. update scene, camera, materials or renderer parameters)
	///

	if (this.animFrame==0)
	{
		this.advanceFrame = false;
	}

	// Advance scene state to next anim frame, if we just exported a rendered frame
	if (this.advanceFrame)
	{	
		gui.sync();
		let no_recompile = true;
		let reinit = false;
		renderer.reset(no_recompile, reinit);
		this.advanceFrame = false;
	}
}

Scene.prototype.postframeCallback = function(snelly, gl)
{
	return;
	let renderer  = snelly.getRenderer();
	let targetSPP = 250.0;

	// User code to post webGL framebuffer data to local server for processing
	if (this.animFrame>=0 && renderer.spp>=targetSPP && this.animFrame<=this.endFrame)
	{
		console.log(this.animFrame);
		let dataURI = gl.canvas.toDataURL();
		var mimetype = dataURI.split(",")[0].split(':')[1].split(';')[0];
		var byteString = atob(dataURI.split(',')[1]);
		var u8a = new Uint8Array(byteString.length);
		for (var i = 0; i < byteString.length; i++) 
		{
			u8a[i] = byteString.charCodeAt(i);
		}
		let blob = new Blob([u8a.buffer], { type: mimetype });
		let r = new XMLHttpRequest();
		r.open('POST', 'http://localhost:3999/' + this.animFrame, false);
		r.send(blob);

		this.advanceFrame = true;
		this.animFrame++;
	}
}
*/

// 2) Launch this server via (assuming node.js is installed):
//     node server.js <output directory for frame images>

// 3) Launch the page to trigger rendering of frames.

// 4) Convert resulting frames to a movie via e.g.:
//    ffmpeg -r 24 -f image2 -i %05d.png -vb 20M video2.mp4

// 5) Optionally, add an audio track via e.g. Avidemux, and export with with options:
//  Video Output:  Mpeg4 AVC (x264)
//  Audio Output:  AAC (Faac)
// Output format:  MP4 Muxer


function getCurrentDirectory() 
{ 
	var fullPath = __dirname; 
	return fullPath; 
}

function ensureDirectoryExistence(dirPath) 
{
	if (fs.existsSync(dirPath)) return true;
	fs.mkdirSync(dirPath);
}

var port = 3999;
var http = require('http');
var fs = require('fs');
var path = require('path');

var outputDirPath = getCurrentDirectory();
if (process.argv.length > 2)
{
	var relDir = process.argv[2];
	outputDirPath = path.join(outputDirPath, relDir);
}

ensureDirectoryExistence(outputDirPath)

http.createServer( function (req, res) 
{
    res.writeHead(200, {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Headers': 'Content-Type, X-Requested-With'
    });
    if (req.method === 'OPTIONS') 
    {
        // Handle OPTIONS requests to work with JQuery and other libs that cause preflighted CORS requests.
        res.end();
        return;
    }
    var idx = req.url.split('/').pop();
    var filename = ("0000" + idx).slice(-5)+".png";
    var filepath = path.join(outputDirPath, filename);

    var img = new Buffer('');
    req.on('data', function(chunk) 
    {
        img = Buffer.concat([img, chunk]);
    });
    req.on('end', function() 
    {
    	console.log('Doing writeFileSync to ' + filepath);
        var f = fs.writeFileSync(filepath, img);
        console.log('Wrote ' + filepath);
        res.end();
    });
}).listen(port, '127.0.0.1');

console.log('Server running at http://127.0.0.1:' + port + '/');
console.log('  files will output to dir ' + outputDirPath)
