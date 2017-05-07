
import os, math
import urllib
import numpy as np

# Script to download refractive index + absorption data from refractiveindex.info
# for a variety of metals, and convert to javascript code

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def frange(x, y, jump):
  while x <= y:
    yield x
    x += jump

def download(metalname, url, Nresample, lmin, lmax):

	response = urllib.urlopen(url)
	table = response.read()
	lines = table.split('\n')
	n = []
	k = []
	for l in lines:
		d = l.split()
		if len(d)>1:
			x = d[0]
			y = d[1]
			if not is_number(x): dataset=y
			else:
				sample = (x, y)
				if   dataset=='n': n.append(sample)
				elif dataset=='k': k.append(sample)
				else: print "don't know dataset '%s'" % dataset; quit()
	 
	if len(n) != len(k):
		print 'bad parsing, mismatched lengths: %d n, %d k' % (len(n), len(k))
	else:
		# resample array into 64 uniformly spaced samples
		w_lst = []
		n_lst = []
		k_lst = []
		for l in range(0, len(n)):
			n_sample = n[l]
			k_sample = k[l]
			if n_sample[0] != k_sample[0]:
				print 'n, k sample wavelength inconsistency! %f %f' % (n_sample[0], k_sample[0]); quit()
			w_lst.append(1000.0*float(n_sample[0]))
			n_lst.append(float(n_sample[1]))
			k_lst.append(float(k_sample[1]))

		dnm = (lmax-lmin)/float(Nresample-1)
		wavelengths_nm = list(frange(lmin, lmax, dnm))
		n_interp = np.interp(wavelengths_nm, w_lst, n_lst)
		k_interp = np.interp(wavelengths_nm, w_lst, k_lst)
		assert len(wavelengths_nm) == Nresample
		assert len(n_interp)       == Nresample
		assert len(k_interp)       == Nresample
		w_str = ",".join(map(str, wavelengths_nm))
		n_str = ",".join(map(str, n_interp))
		k_str = ",".join(map(str, k_interp))
	
		# wavelength_nm: new Float32Array([%s]),
		code = '''function tabulated_%s() { // %s
	return { n: new Float32Array([%s]),
			 k: new Float32Array([%s]) } }
''' % (metalname, "%s samples of n, k between %fnm and %fnm"%(Nresample, lmin, lmax), n_str, k_str)
		
		return code
	

print download("copper",     "https://refractiveindex.info/tmp/main/Cu/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("silver",     "https://refractiveindex.info/tmp/main/Ag/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("gold",       "https://refractiveindex.info/tmp/main/Au/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("aluminium",  "https://refractiveindex.info/tmp/main/Al/Rakic.txt",   Nresample=64, lmin=390.0, lmax=750.0)
print download("chromium",   "https://refractiveindex.info/tmp/main/Cr/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("iron",       "https://refractiveindex.info/tmp/main/Fe/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("molybdenum", "https://refractiveindex.info/tmp/main/Mo/Ordal.txt",   Nresample=64, lmin=390.0, lmax=750.0)
print download("nickel",     "https://refractiveindex.info/tmp/main/Ni/Ordal.txt",   Nresample=64, lmin=390.0, lmax=750.0)
print download("lead",       "https://refractiveindex.info/tmp/main/Pb/Ordal.txt",   Nresample=64, lmin=390.0, lmax=750.0)
print download("palladium",  "https://refractiveindex.info/tmp/main/Pd/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("platinum",   "https://refractiveindex.info/tmp/main/Pt/Werner.txt",  Nresample=64, lmin=390.0, lmax=750.0)
print download("silicon",    "https://refractiveindex.info/tmp/main/Si/Aspnes.txt",  Nresample=64, lmin=390.0, lmax=750.0)
print download("titanium",   "https://refractiveindex.info/tmp/main/Ti/Johnson.txt", Nresample=64, lmin=390.0, lmax=750.0)
print download("tungsten",   "https://refractiveindex.info/tmp/main/W/Ordal.txt",    Nresample=64, lmin=390.0, lmax=750.0)
print download("zinc",       "https://refractiveindex.info/tmp/main/Zn/Werner.txt",  Nresample=64, lmin=390.0, lmax=750.0)


