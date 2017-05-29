

// Launch via (assuming node.js is installed):
//  node server.js <output directory for images>


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
