
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
    const INF = 1e20;
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

window.MathJax = {

    jax: ["input/TeX", "output/SVG"],
    extensions: ["tex2jax.js", "MathMenu.js", "MathZoom.js"],
    showMathMenu: false,
    showProcessingMessages: false,
    messageStyle: "none",
    SVG: {
        useGlobalCache: false
    },
    TeX: {
        extensions: ["AMSmath.js", "AMSsymbols.js", "autoload-all.js"]
    },

    AuthorInit: function() {
        MathJax.Hub.Register.StartupHook("End", function() {
        var mj2img = function(texstring, callback) {
            var input = texstring;
            var wrapper = document.createElement("div");
            wrapper.innerHTML = input;
            var output = { svg: "", img: ""};
            MathJax.Hub.Queue(["Typeset", MathJax.Hub, wrapper]);
            MathJax.Hub.Queue(function() {
            var mjOut = wrapper.getElementsByTagName("svg")[0];
            mjOut.setAttribute("xmlns", "http://www.w3.org/2000/svg");
            // thanks, https://spin.atomicobject.com/2014/01/21/convert-svg-to-png/
            output.svg = mjOut.outerHTML;
            output.svg.width = 512;
            output.svg.height = 128;
            var image = new Image();
            image.width = SDF_MATH_WIDTH;
            image.height = SDF_MATH_HEIGHT;
            image.src = 'data:image/svg+xml;base64,' + window.btoa(unescape(encodeURIComponent(output.svg)));
            image.onload = function() {
                callback(image);
            };
            });
        }

        mj2img(SDF_MATH_TEXT, function(image) {
            
            var img = image;
            this.W = img.width;
            this.H = img.height;

            // temporary arrays for the distance transform
            this.gridOuter = new Float64Array(this.W * this.H);
            this.gridInner = new Float64Array(this.W * this.H);
            this.f = new Float64Array(Math.max(this.W, this.H));
            this.d = new Float64Array(Math.max(this.W, this.H));
            this.z = new Float64Array(Math.max(this.W, this.H) + 1);
            this.v = new Int16Array(Math.max(this.W, this.H));

            this.canvas = document.createElement('canvas');
            this.canvas.width = this.W;
            this.canvas.height = this.H;
            this.ctx = this.canvas.getContext('2d');

            const INF = 1e20;
            this.ctx.drawImage(img, 0, 0, img.width, img.height);
            var imgData = this.ctx.getImageData(0, 0, this.W, this.H);
            for (var i = 0; i < this.W * this.H; i++) {
                var a = imgData.data[i * 4 + 3] / 255; // alpha value
                this.gridOuter[i] = a === 1 ? 0 : a === 0 ? INF : Math.pow(Math.max(0, 0.5 - a), 2);
                this.gridInner[i] = a === 1 ? INF : a === 0 ? 0 : Math.pow(Math.max(0, a - 0.5), 2);
            }
            edt(this.gridOuter, this.W, this.H, this.f, this.d, this.v, this.z);
            edt(this.gridInner, this.W, this.H, this.f, this.d, this.v, this.z);

            var sdf = new Float32Array(this.W * this.H);
            for (i = 0; i < this.W * this.H; i++) {
                var d = this.gridOuter[i] - this.gridInner[i];
                sdf[i] = d / this.W; // SDF in units of texture width
            }
            MATHSDF_GENERATED_CALLBACK(sdf, this.W, this.H);

        });

    });
  }
};

