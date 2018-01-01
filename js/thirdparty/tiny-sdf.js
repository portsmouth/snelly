const INF = 1e20;

function TextSDF(fontSize, resHorizontal, resVertical, fontFamily, fontWeight) {
    this.fontSize = fontSize;
    this.fontFamily = fontFamily || 'serif';
    this.fontWeight = fontWeight || 'normal';

    this.W = resHorizontal;
    this.H = resVertical;

    this.canvas = document.createElement('canvas');
    this.canvas.width = this.W;
    this.canvas.height = this.H;

    this.ctx = this.canvas.getContext('2d');
    this.ctx.font = this.fontWeight + ' ' + this.fontSize + 'px ' + this.fontFamily;
    this.ctx.textBaseline = 'middle';
    this.ctx.fillStyle = 'black';

    // temporary arrays for the distance transform
    this.gridOuter = new Float64Array(this.W * this.H);
    this.gridInner = new Float64Array(this.W * this.H);
    this.f = new Float64Array(Math.max(this.W, this.H));
    this.d = new Float64Array(Math.max(this.W, this.H));
    this.z = new Float64Array(Math.max(this.W, this.H) + 1);
    this.v = new Int16Array(Math.max(this.W, this.H));
}

TextSDF.prototype.getResHorizontal = function() 
{
    return this.W;
}

TextSDF.prototype.getResVertical = function() 
{
    return this.H;
}

TextSDF.prototype.draw = function (char) {
    this.ctx.clearRect(0, 0, this.W, this.H);
    this.ctx.fillText(char, 0, this.H/2, this.W);

    var imgData = this.ctx.getImageData(0, 0, this.W, this.H);
    var alphaChannel = new Float32Array(this.W * this.H);

    for (var i = 0; i < this.W * this.H; i++) {
        var a = imgData.data[i * 4 + 3] / 255; // alpha value
        this.gridOuter[i] = a === 1 ? 0 : a === 0 ? INF : Math.pow(Math.max(0, 0.5 - a), 2);
        this.gridInner[i] = a === 1 ? INF : a === 0 ? 0 : Math.pow(Math.max(0, a - 0.5), 2);
    }

    edt(this.gridOuter, this.W, this.H, this.f, this.d, this.v, this.z);
    edt(this.gridInner, this.W, this.H, this.f, this.d, this.v, this.z);

    for (i = 0; i < this.W * this.H; i++) {
        var d = this.gridOuter[i] - this.gridInner[i];
        alphaChannel[i] = d / this.W; // SDF in units of texture width
    }

    return alphaChannel;
};

// 2D Euclidean distance transform by Felzenszwalb & Huttenlocher https://cs.brown.edu/~pff/dt/
function edt(data, width, height, f, d, v, z) {
    for (var x = 0; x < width; x++) {
        for (var y = 0; y < height; y++) {
            f[y] = data[y * width + x];
        }
        edt1d(f, d, v, z, height);
        for (y = 0; y < height; y++) {
            data[y * width + x] = d[y];
        }
    }
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            f[x] = data[y * width + x];
        }
        edt1d(f, d, v, z, width);
        for (x = 0; x < width; x++) {
            data[y * width + x] = Math.sqrt(d[x]);
        }
    }
}

// 1D squared distance transform
function edt1d(f, d, v, z, n) {
    v[0] = 0;
    z[0] = -INF;
    z[1] = +INF;

    for (var q = 1, k = 0; q < n; q++) {
        var s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        while (s <= z[k]) {
            k--;
            s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = +INF;
    }

    for (q = 0, k = 0; q < n; q++) {
        while (z[k + 1] < q) k++;
        d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
    }
}
