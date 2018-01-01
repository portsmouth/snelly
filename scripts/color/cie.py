#!/usr/bin/python
import csv
from decimal import Decimal
import numpy

#############################################################################
# Script to generate wavelengthToXYZ.js
#############################################################################

def rgbToXyz(RGB):
    XYZ = [0, 0, 0]
    XYZ[0] = 0.4124564*RGB[0] + 0.3575761*RGB[1] + 0.1804375*RGB[2]
    XYZ[1] = 0.2126729*RGB[0] + 0.7151522*RGB[1] + 0.0721750*RGB[2]
    XYZ[2] = 0.0193339*RGB[0] + 0.1191920*RGB[1] + 0.9503041*RGB[2]
    return XYZ

def xyz_to_spectrum(XYZ):
    c = [0, 0, 0]
    c[0] =  3.38214566*XYZ[0] - 2.58540997*XYZ[1] - 0.40649004*XYZ[2]
    c[1] = -2.58540997*XYZ[0] + 3.20943158*XYZ[1] + 0.22767094*XYZ[2]
    c[2] = -0.40649004*XYZ[0] + 0.22767094*XYZ[1] + 0.70334476*XYZ[2]
    return c

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
    X.append( float(row[1]) )
    Y.append( float(row[2] ) )
    Z.append( float(row[3] ) )


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
# Also compute cross-term matrix needed to map texture colors into XYZ basis function coefficients
c_x = 0.0;
c_y = 0.0;
c_z = 0.0;
c_xx = 0.0; c_xy = 0.0; c_xz = 0.0
c_yx = 0.0; c_yy = 0.0; c_yz = 0.0
c_zx = 0.0; c_zy = 0.0; c_zz = 0.0
data = ''
for n in range(0, N):

    l = l_resampled[n]
    x = float(X_resampled[n])
    y = float(Y_resampled[n])
    z = float(Z_resampled[n])

    c_x  +=   x * dl
    c_y  +=   y * dl
    c_z  +=   z * dl
    c_xx += x*x * dl
    c_xy += x*y * dl
    c_xz += x*z * dl
    c_yx += y*x * dl
    c_yy += y*y * dl
    c_yz += y*z * dl
    c_zx += z*x * dl
    c_zy += z*y * dl
    c_zz += z*z * dl

    data += '%f, %f, %f, 0.0' % (x, y, z)
    if n!=N-1: data += ', '
    if (n+1)%4==0: data += '\n\t\t'

ns = [0, 512, 1023]
for n in ns:
    print '%d: %f %f %f' % (n, float(X_resampled[n]), float(Y_resampled[n]), float(Z_resampled[n]))

#print '''
#// A table of %d vec4 tristimulus XYZ values, corresponding to the %d wavelength samples
#// between 390.0 and 750.0 nanometres
#function wavelengthToXYZTable() {
#    return new Float32Array([
#            %s
#        ]);
#}
#''' % (N, N, data)

cm_norm = c_x/(750.0 - 390.0)
print 'XYZ norms: \n', cm_norm, c_y/(750.0 - 390.0), c_z/(750.0 - 390.0)

M = numpy.matrix( [[c_xx, c_xy, c_xz],
                          [c_yx, c_yy, c_yz],
                          [c_zx, c_zy, c_zz]])
print '\nXYZ correlation matrix: \n', M
print '\nXYZ correlation matrix inverse: \n', M.I

# scale by an ad-hoc factor to make max RGB component about 1
xyzToSpc = M.I
print '\nXYZ to SPC'
print xyzToSpc

rgbToXyz = numpy.matrix( [[0.4124564, 0.3575761, 0.1804375],
                          [0.2126729, 0.7151522, 0.0721750],
                          [0.0193339, 0.1191920, 0.9503041]])
print '\nRGB to SPC'
rgbToSpc = xyzToSpc * rgbToXyz
white_unscaled = rgbToSpc * numpy.matrix([[1.0], [1.0], [1.0]])
scale = (1.0/cm_norm) / numpy.sum(white_unscaled)
rgbToSpc = scale * rgbToSpc
print rgbToSpc

#xyz_to_spectrum = numpy.matrix( [[3.38214566, -2.58540997, -0.40649004],
#                                 [-2.58540997, 3.20943158,  0.22767094],
#                                 [-0.40649004, 0.22767094,  0.70334476]])
#print xyz_to_spectrum * rgbToXyz
