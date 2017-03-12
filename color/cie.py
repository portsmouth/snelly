
import csv
from decimal import Decimal
import numpy

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

# Convert to RGB values using sRGB transformation:
data = ''
for n in range(0, N):

	l = l_resampled[n]
	x = float(X_resampled[n])
	y = float(Y_resampled[n])
	z = float(Z_resampled[n])

	r =  3.2406*x - 1.5372*y - 0.4986*z
	g = -0.9689*x + 1.8758*y + 0.0415*z
	b =  0.0557*x - 0.2040*y + 1.0570*z

	data += '%f, %f, %f, 0.0' % (r, g, b)
	if n!=N-1: data += ', ' 
	if (n+1)%4==0: data += '\n\t\t'

print '''
// A table of %d vec4 sRGB colors, corresponding to the %d wavelength samples
// between 390.0 and 750.0 nanometres
function wavelengthToRgbTable() {
    return new Float32Array([
    		%s
        ]);
}
''' % (N, N, data)