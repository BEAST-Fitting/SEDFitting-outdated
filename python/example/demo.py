""" DEMO Script, aiming at presented a quick tutorial rather than performance
This script verify that all necessary files are available
Python package imports are done when needed for reading clarity (not in the pure
pythonic philosophy)
"""

""" 
============================ DEMO INTEGRITY  ===================================
For the purpose of this script, we want to be sure you have downloaded the required
libraries
(this function is called first during the execution)
"""
import os, inspect
import anased
demo_data_file = './example.obs.csv'
def check_libs_requirements():
	""" Check if you have all the files to run this demo 
	"""
	root = '/'.join(os.path.abspath(inspect.getfile(anased)).split('/')[:-1])

	libs = {'filters': 'filters.hd5',
		'sedgrid': 'PHATSEDs_basel_padovaiso.fits',
		'vega': 'vega.hd5'}

	for k in libs:
		assert (os.path.isfile(root+'/libs/'+libs[k])), 'Missing %s (%s).\n   run "> make libs"' % (k.upper(), libs[k])
	
	
	'Generating DEMO DATA FILE (%s).' % (demo_data_file)
	generateFakeData(output=demo_data_file, N=10, Av_max=0.1)
                
""" 
================================ INPUT ====================================
define inputs from a "standard" table, i.e, readable by
anased.tools.mytables.load(), pretty much any common format (csv, tsv, ascii,
fits, hdf5)
"""
# for this example, we derive the provided Observation class.
# For this particular example, the dataset is given as a csv file, assuming:
#    * measurements in 6 filters
#    * 6 magnitude measurements, in VEGA mag system.
#    * 6 associated uncertainties (Gaussian! for the likelihood)
#    * Possibly many other columns, e.g. RA, Dec, ...
#    * bad values are set to any value above 99. mag or nan or inf

#normalized names, basically the norm is <FACILITY>_<INSTRUMENT>_<FILTERNAME> (in caps) 
filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',
		'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']   


#Data are in Vega magnitudes, we will do conversions on the fly, but we need
#some info about Vega first (reason why with require VEGA (vega.hd5) library)
from anased.tools.vega import Vega, from_Vegamag_to_Flux
with Vega() as v:
	vega_f, vega_mag, lamb = v.getMag(filters)

#Define a standard observation interface to our data file called PhatData
from anased.observations import Observations
class PhatData(Observations):
   """ PHAT catalog for stars in M31

	    The original Observation class already handles:
	    	* bad values
		* distance modulus / distance
		* and load the data from the file
   """
   def __init__(self, inputFile, distanceModulus=0.):
	""" Class Contructor 
		inputFile	str	filename incl. path
		distanceModulus float	apply a distance modulus (default 24.3mag)
	"""
	desc = 'PHAT catalog: %s' % inputFile
	Observations.__init__(self, inputFile,distanceModulus, desc=desc)
	self.setFilters( filters )
	self.setBadValue(99.0)

   @from_Vegamag_to_Flux(lamb, vega_mag)
   def getObs(self, num):
	""" Using the decorator @from_Vegamag_to_Flux 
	    to return values in flux (not in flux/flux_vega) from vega
	    magnitudes

	    num		int	line in the table (0..n)
	
	    By default, Observation.getObs() will consider:
	    	* 6 columns with names corresponding to the filters
		* 6 columns with the filternames+'err' as the measurement errors
	    
	    Bad values are already handled by default.

	    Returns the fluxes, errors and mask of observation 'num'.
	    Order matters!! (filters order)
	"""
	return Observations.getObs(self, num)

   # you can define any other function you want as long as it was not already
   # defined in Observations
   def getObsinMag(self, num):
	""" Returns the original catalog magnitudes
		Looks like the same as getObs, but not decorated by
   			@from_Vegamag_to_Flux(lamb, vega_mag)
		Which will not convert the values to fluxes
	"""
	return Observations.getObs(self, num)


""" 
============================ MODEL GRID ===================================
Load the models.
We assume the models are SEDs, defined for the same filters (in the same
order)

Each models has a set of 
   * 6 absolute fluxes 
	in units that corresponds to observed fluxes, i.e., in erg/s/cm2/AA
	(We convert vega mags to the same fluxes units)
   * many internal properties
   	e.g., age, mass, temperature, gravity, luminosity ...

We will use a prepared SED grid (assuming someone did correctly the preparation)
that follows the standards of the code.
"""
# I do the two lines below in the main function instead
#    Because you may want to consider the model grid as a variable source
#from anased import grid
#sedGrid = grid.FileSpectralGrid('libs/SEDs_basel_padovaiso.fits')


"""
============================ MAIN CODE ===================================
"""

""" Define the actual Job """
from anased import getFluxAttenuation, computeLogLikelihood
from numpy import exp
# getFluxAttenuation check the variables and the units (shortcut)
def task(lamb, f, e, m, fm, extLaw, **kwargs):
	""" Shortcut to compute the log likelihood of the SED with the models
	    for a given extinction parameter set.
	INPUTS:
		lamb	np.ndarray[float, ndim=1]	array of wavelengths in AA
		f	np.ndarray[float, ndim=1]	array of fluxes
		e	np.ndarray[float, ndim=1]	array of flux errors
		m	np.ndarray[bool, ndim=1]	mask array to apply during the calculations
						        mask.shape = flux.shape
		fm	np.ndarray[float, ndim=2]	array of modeled fluxes (Nfilters , Nmodels)
		extLaw	extinction.ExtinctionLaw 	instance of extinction law
	KEYWORDS:
		**kwargs is forwarded to the getFluxAttenuation call
	OUTPUTS:
		ln(L) 	np.ndarray[float, ndim=1]	array of ln(L) values (Nmodels)

	Note: this function is a copy of anased.job for clarity of the computations
	"""
	# get attetuation values
	#tau = extLaw.function( lamb * 1e-4, Alambda = False, **kwargs)
	tau = getFluxAttenuation(extLaw, lamb, **kwargs)
	#deredden the observed flux (faster than adding reddening to all models
	deredflux = f*exp(tau)
	#compute lnp
	lnp = computeLogLikelihood(deredflux, e, fm, normed=False, mask=m)	
	return lnp

import numpy
""" Do some figures """
from matplotlib.ticker import MaxNLocator
import pylab as plt
def plotPDF(r, sedGrid, Av, Q = 'logg logT logL logM logA Av Z', starnum=1):
	"""
	Plot the marginal pdfs of some quantities
	INPUTS:
		r	ndarray[float, ndim=2]	lnp array
		sedGrid grid.ModelGrid		grid of models associated to lnp
		Av	ndarray[float, ndim=1]	list of Av values
	
	KEYWORDS:
		Q	string		list of quantities to plot (separated by spaces)
		starnum str/int		number used in the figure filename
	
	plots are saved into png files in the current directory.
	"""
	_q = Q.split()
	_r = numpy.exp(r)
	_r /= _r.sum()

	d = {}
	r1 = _r.sum(1)
	rAv = _r.sum(0)

	N  = len(_q)
	nc = 3 
	nl = (N - N//nc)

	plt.figure(figsize=(10,10))
	j=0
	for k in _q:
		ax = plt.subplot(nl, nc, j+1)
		if k.lower() != 'av': 
			n, b = numpy.histogram(sedGrid.grid[k], bins=20, weights=r1)
		else:
			n, b = numpy.histogram(Av, bins=len(Av), weights=rAv)
		ax.step(b[:-1], n.astype(float)/n.sum(), color='black', lw=2, where='post')
		plt.xlabel(k)
		plt.ylabel('P(data |' + k+')')
		ax.xaxis.set_major_locator(MaxNLocator(5))
		ax.yaxis.set_major_locator(MaxNLocator(4))
		j+=1

	plt.subplots_adjust(hspace=0.4, wspace=0.4)
	plt.savefig('stellarpdf_'+str(starnum)+'.png', format='png', dpi=150, bbox_inches='tight')
	plt.close()


def generateFakeData(filters=filters, output='example.obs.csv', N=20, Av_min=0., Av_max=3., Rv=3.1, err=0.05):
	""" Function that generates fake data from the grid of models """
	from anased import grid
	from anased.extinction import Cardelli
	from anased import getFluxAttenuation
	from anased.tools.vega import Vega
	import numpy

	with Vega() as v:
		vega_f, vega_mag, lamb = v.getMag(filters)

	g      = grid.FileSpectralGrid('anased/libs/PHATSEDs_basel_padovaiso.fits')
	oAv    = Cardelli()
	lamb   = g.lamb

	fakein = numpy.random.randint(0, g.grid.nrows, N)
	Av0    = numpy.random.uniform(Av_min, Av_max, N)
	magerr = numpy.array( [0.05]*len(filters) )

	d = g.grid[fakein]
	d.addCol(Av0, name='Av0')
	d.addCol(fakein, name='modelId')
	fakemag = numpy.empty((N, len(filters)),dtype=float)
	for k in range(N):
		fakesed       = numpy.copy(g.seds[fakein[k],:])
		tau           = getFluxAttenuation(oAv, lamb, Av = Av0[k], Rv = Rv)
		fakesed      *= exp(-tau)
		fakemag[k,:]  = -2.5*numpy.log10(fakesed) - vega_mag
			
	for k in range(len(filters)):
		d.addCol( fakemag[:,k] , name=filters[k] )
		d.addCol( err*numpy.ones(N,dtype=float), name=filters[k]+'err' )

	d.header['NAME'] = 'Fake data'
	d.header['COMMENT'] = 'values are in Vega Magnitudes'
	d.write(output)	



def main():
	obs = PhatData(demo_data_file, distanceModulus = 0.)

	# we suppose a Cardelli extinction law
	from anased.extinction import Cardelli
	av_law  = Cardelli()
	Av_min  = 0.
	Av_max  = 3.
	Av_step = 0.3 # faster :p
	Rv      = 3.1 # fixed in this example
	Av      = numpy.arange(Av_min,Av_max+Av_step, Av_step)

	#load the models
	from anased import grid
	sedGrid = grid.FileSpectralGrid('anased/libs/PHATSEDs_basel_padovaiso.fits')
	lamb    = sedGrid.lamb

	#define where to keep the results
	for idx in range(len(obs)):
		print 'Taking care of obs. %d / %d' % (idx, len(obs))
		f, e, m = obs.getObs(idx)
		r = numpy.empty( (sedGrid.seds.shape[0], len(Av)), dtype=float )
		for k in range(len(Av)):
			r[:, k] = task(lamb[:], f, e, m, numpy.copy(sedGrid.seds),
					av_law, Av=Av[k], Rv=3.1) 	
		plotPDF(r, sedGrid, Av, starnum=idx)

if __name__ == '__main__':
	check_libs_requirements()
	main()

