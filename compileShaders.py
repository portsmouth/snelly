
import os
import imp

###################################
# compile shader source
###################################

ShaderDir = 'shaders'
src = "var Shaders = {\n\n"

commonPath = ShaderDir + '/common.glsl'
commonCode = open(commonPath).read().strip().split('\n')

debug = False

for f in os.listdir(ShaderDir):

    if f == 'common.glsl': continue;
    if debug: 
        print '\n' + '*'* 80
        print f
        print '*'* 80 + '\n'
    if f.find('.glsl') == -1: continue

    name = f.replace('.glsl', '')
    src += "'%s': `\n" % name

    path = os.path.join(ShaderDir, f)
    lines = open(path).read().strip().split('\n')
    
    code = list(commonCode)
    code.append('')
    code.extend(lines)

    for i, line in enumerate(code):
        src += line + '\n'
        if debug: print line
        if i==len(code)-1:
            src += "`,\n\n"

src += "}"

open("js/shaders.js", 'w').write(src)

