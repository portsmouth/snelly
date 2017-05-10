
import csv
from decimal import Decimal
import numpy

#############################################################################
# Script to generate wavelengthToXYZ.js
#############################################################################

# 2-deg XYZ tristimulus CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals
# (from http://www.cvrl.org/cmfs.htm)
csvfile = open("lin2012xyz2e_fine_7sf.csv", 'rb')
cie = csv.reader(csvfile, delimiter=',')

l_data = []
X = []
Y = []
Z = []
for row in cie:
	l_data.append( Decimal(row[0]) )
	X.append( Decimal(row[1]) )
	Y.append( Decimal(row[2]) )
	Z.append( Decimal(row[3]) )


# Convert to 1024 XYZ values sampled at wavelengths between 390.0 and 750.0 nm
N = 1024
l_resampled = []
dl = (750.0 - 390.0)/float(N-1)
for n in range(0, N):
	l_resampled.append(390.0 + dl*n)

X_resampled = numpy.interp(l_resampled, l_data, X)
Y_resampled = numpy.interp(l_resampled, l_data, Y)
Z_resampled = numpy.interp(l_resampled, l_data, Z)

# Convert to RGB values using sRGB transformation.
# Also compute the cross-term matrix needed to map sRGB texture colors into RGB basis function coefficients
c_rr = 0.0; c_rg = 0.0; c_rb = 0.0
c_gr = 0.0; c_gg = 0.0; c_gb = 0.0
c_br = 0.0; c_bg = 0.0; c_bb = 0.0
data = ''
for n in range(0, N):

	l = l_resampled[n]
	x = float(X_resampled[n])
	y = float(Y_resampled[n])
	z = float(Z_resampled[n])

	# XYZ to RGB transformation
	r =  3.2406*x - 1.5372*y - 0.4986*z
	g = -0.9689*x + 1.8758*y + 0.0415*z
	b =  0.0557*x - 0.2040*y + 1.0570*z

	c_rr += r*r * dl
	c_rg += r*g * dl
	c_rb += r*b * dl
	c_gr += g*r * dl
	c_gg += g*g * dl
	c_gb += g*b * dl
	c_br += b*r * dl
	c_bg += b*g * dl
	c_bb += b*b * dl

	data += '%f, %f, %f, 0.0' % (r, g, b)
	if n!=N-1: data += ', ' 
	if (n+1)%4==0: data += '\n\t\t'

print '''
// A table of %d vec4 RGB values in sRGB color space, corresponding to the %d wavelength samples
// between 390.0 and 750.0 nanometres
function wavelengthToRGBTable() {
    return new Float32Array([
    		%s
        ]);
}
''' % (N, N, data)


M = numpy.matrix( [[c_rr, c_rg, c_rb], 
	               [c_gr, c_gg, c_gb], 
	               [c_br, c_bg, c_bb]]) 
print M.I


