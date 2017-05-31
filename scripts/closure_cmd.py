
import os, io
import rjsmin
 
# Construct the command for Google Closure Compiler, to generate minified code

jsDestPath = "../js/compiled/snelly.min.js"

jsSources = [ 
'../js/thirdparty/jquery-1.11.3.min.js',
'../js/thirdparty/three/three.min.js',
'../js/thirdparty/three/libs/stats.min.js',
'../js/thirdparty/three/libs/dat.gui.min.js',
'../js/thirdparty/three/controls/OrbitControls.js',
'../js/gl.js',
'../js/gui.js',
'../js/shaders.js',
'../js/color.js',
'../js/materials.js',
'../js/spectra.js',
'../js/pathtracer.js',
'../js/snelly.js' 
]

files = ' '.join(jsSources)

cmd = '''
java -jar closure-compiler/compiler.jar --js %s -O WHITESPACE_ONLY --js_output_file %s
''' % (files, jsDestPath)

print cmd