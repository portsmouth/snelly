
import os

ShaderDir = 'shaders'
src = "var Shaders = {\n\n"

for file in os.listdir(ShaderDir):

    if file.find('.glsl') == -1:
        continue

    name = file.replace('.glsl', '')
    src += "'%s': \n" % name

    path = os.path.join(ShaderDir, file)
    lines = open(path).read().strip().split('\n')

    for i, line in enumerate(lines):
        src += "\t'" + line + "\\n'"
        if i<len(lines)-1:
            src += " +\n"
        else:
            src += ",\n\n"

src += "}"

open("js/shaders.js", 'w').write(src)