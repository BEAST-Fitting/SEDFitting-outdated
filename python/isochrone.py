"""
Isochrone class

Intent to implement a generic module to manage isochrone mining from various
sources.


"""
import warnings
import numpy
from numpy import interp
from numpy import log10
from tools import units
import inspect, os
import tables
from tools import mytables

localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])


class Isochrone(object):
	def __init__(self, name='', *args, **kwargs):
		self.name = name
	def metalToFeH(self,metal):
		""" Convert Z to [Fe/H] values 
			Zsun = 0.02 <-> [Fe/H]sun = -4.33
			Z = [ 0.0004, 0.004, 0.008, 0.02, 0.05 ] 
		   [Fe/H] = [ -1.7  , -0.7 , -0.4 , 0   , 0.4  ]
		"""
		return numpy.log10(metal/0.02)

	def FeHtometal(self,feh):
		""" Convert Z to [Fe/H] values 
			Zsun = 0.02 <-> [Fe/H]sun = -4.33
			Z = [ 0.0004, 0.004, 0.008, 0.02, 0.05 ] 
		   [Fe/H] = [ -1.7  , -0.7 , -0.4 , 0   , 0.4  ]
		"""
		return 10**feh * 0.02
	
	def _get_isochrone(self, *args, **kwargs):
		""" Retrieve isochrone from the original source 
			internal use to adapt any library
		"""
		pass

	def _get_continuous_isochrone(self, *args, **kwargs):
		""" Return a resampled isochrone accounting for variations
			useful for continuous sampling
		"""
		dm = kwargs.pop('dm', 0.01)
		dt = kwargs.pop('dt', 0.01)
		dl = kwargs.pop('dl', 0.01)

		iso = self._get_isochrone(*args, **kwargs)
		logT, logg, logL, logM = iso['logT'], iso['logg'], iso['logL'], iso['logM']

		#compute vector of discrete derivaties for each quantity
		dlogm = numpy.abs(numpy.diff((logM)))
		dlogT = numpy.abs(numpy.diff(logT))
		dlogg = numpy.abs(numpy.diff(logg))
		#set up vectors for storage
		newm= []
		newt = []
		newg = []
		newl = []

		#define the maximum allowable difference between points
		maxdiff = 0.01
		for i in range(len(dlogm)):
			#check to see if difference for point i has a problem in some dimension
			if ((dlogm[i] >= dm)| (dlogT[i] >= dt) | (dlogg[i] >= dl)):
				#compute the number of points
				npts = numpy.int(dlogm[i]/maxdiff) + numpy.int(dlogT[i]/maxdiff) + numpy.int(dlogg[i]/maxdiff)
				#construct new 1d grids in each dimension, being careful about endpoints
				#append them to storage vectors
				newm = numpy.append(newm, numpy.linspace(logM[i], logM[i+1], npts+1, endpoint=False))
				newt = numpy.append(newt, numpy.linspace(logT[i], logT[i+1], npts+1, endpoint=False))
				newg = numpy.append(newg, numpy.linspace(logg[i], logg[i+1], npts+1, endpoint=False))
				newl = numpy.append(newl, numpy.linspace(logL[i], logL[i+1], npts+1, endpoint=False))
			else: #if the maximumum allowable difference is small, then just store the good point
				newm = numpy.append(newm, logM[i])
				newt = numpy.append(newt, logT[i])
				newg = numpy.append(newg, logg[i])
				newl = numpy.append(newl, logL[i])
		#tack on the last point on the grid, as the loop is one element short
		newm = numpy.append(newm, logM[-1])
		newt = numpy.append(newt, logT[-1])
		newg = numpy.append(newg, logg[-1])
		newl = numpy.append(newl, logL[-1])

		data = dict(logM=newm, logT=newt, logg=newg, logL=newl)
		table = mytables.Table(data)

		for k in iso.header.keys():
			table.header[k] = iso.header[k]

		table.header['NAME'] = 'Resampled '+ table.header['NAME']

		table.header['dlogT'] = dt
		table.header['dlogM'] = dm
		table.header['dlogg'] = dl
		
		return table

class padova2010(Isochrone):
        def __init__(self):
		self.name = 'Padova 2010 (Marigo 2008 + Girardi 2010)'
		self.source = localpath+'/libs/padova2010.iso.fits'
		self._load_table_(self.source)
		self.ages = 10**numpy.unique(self.data['logA']) 
		self.Z    = numpy.asarray([0.019]) 

	def _load_table_(self, source):
		t = mytables.load(self.source)
		data = {}
		for k in t.keys():
			data[k] = t[k]
		#Alias columns
		data['logM'] = log10(numpy.asarray(data['M_ini']))
		data['logg'] = numpy.asarray(data['logG'])
		data['logT'] = numpy.asarray(data['logTe'])
		data['logL'] = numpy.asarray(data['logL/Lo'])
		data['logA'] = numpy.asarray(data['log(age/yr)'])
		#clean columns
		data.pop('log(age/yr)')
		data.pop('M_ini')
		data.pop('logG')
		data.pop('logTe')
		data.pop('logL/Lo')

		self.data = mytables.Table(data, name='Isochrone from %s' % self.name)


	def _get_isochrone(self, age, metal=None, FeH=None, inputUnit=units.yr, masses = None, *args, **kwargs):
		""" Retrieve isochrone from the original source 
			internal use to adapt any library
		"""
		if isinstance(age, units.Quantity):
			_age = int(age.rescale('yr'))
		else:
			_age = int(age * inputUnit.rescale('yr'))
		
		assert ((metal != None) | (FeH != None)), "Need a chemical par. value."
		
		if (metal != None) & (FeH != None):
			warnings.warn('both Z & [Fe/H] provided, ignoring [Fe/H].',Warning)

		if metal == None:
			metal = self.FeHtometal(FeH)

		if metal != None:
			warnings.warn('padova library is under construction only Zsun available now.',Warning)

		data = {}
		if _age in self.ages:
			#no interpolation, isochrone already in the file
			t = self.data.selectWhere( 'logA == _age', condvars={'_age':log10(_age)} )
			for kn in t.keys():
				data [kn] = numpy.asarray(t[kn]) 
		else:
			#interpolate between isochrones
			d      = (self.ages - float(_age))**2
			a1, a2 = 10**self.ages[numpy.argsort(d)[:2]]
			#print "Warning: Interpolation between %d and %d Myr" % (a1,a2)
			r = numpy.log10(_age/a1)/numpy.log10(a2/a1)

			t1 =self.data.selectWhere( 'logA == _age', condvars={'_age':log10(a1)} )
			t2 =self.data.selectWhere( 'logA == _age', condvars={'_age':log10(a2)} )

			stop = min(t1.nrows,t2.nrows)
		
			for kn in t1.colnames:
				y2 = t2[kn][:stop]
				y1 = t1[kn][:stop]
				data[kn] = y2 * r + y1 * (1.-r)
				del y1, y2
		

		#mass selection
		if masses != None:
			#masses are expected in logM for interpolation
			if masses.max() > 2.3:
				_m = numpy.log10(masses)
			else:
				_m = masses
			data_logM = data['logM'][:]
			for kn in data:
				data[kn] = interp(_m, data_logM, data[kn])
				
		table = mytables.Table(data, name='Isochrone from %s' % self.name)
		table.header['metal'] = metal
		table.header['time'] = _age
		return table	
		



class pegase(Isochrone):
        def __init__(self):
		self.name   = 'Pegase.2 (Fioc+1997)'
		self.source = localpath+'/libs/pegase.iso.hd5'
		self.data   = tables.openFile(self.source)
		self.ages   = numpy.asarray([k.attrs.time for k in self.data.root.Z02])*1e6
		self.Z      = numpy.asarray([ float('0.'+k[1:]) for k in self.data.root._g_listGroup(self.data.getNode('/'))[0]])

	def __getstate__(self):
		self.data.close()
		self.data = None
		return self.__dict__

	def __setstate__(self, d):
		self.__dict__ = d
		self.data = tables.openFile(self.source)

	def __del__(self):
		if self.data != None:
			self.data.close()

	def _get_isochrone(self, age, metal=None, FeH=None, inputUnit=units.yr, masses = None, *args, **kwargs):
		""" Retrieve isochrone from the original source 
			internal use to adapt any library
		"""
		if isinstance(age, units.Quantity):
			_age = int(age.rescale('Myr'))
		else:
			_age = int(age * inputUnit.rescale('Myr'))
		
		assert ((metal != None) | (FeH != None)), "Need a chemical par. value."
		
		if (metal != None) & (FeH != None):
			warnings.warn('both Z & [Fe/H] provided, ignoring [Fe/H].',Warning)

		if metal == None:
			metal = self.FeHtometal(FeH)

		groups = self.data.root._g_listGroup(self.data.root)[0]
		assert ('Z'+str(metal)[2:] in groups), "Metal %f not find in %s" % (metal, str(groups))
		node = self.data.getNode('/Z'+str(metal)[2:])
		groups = self.data.root._g_listGroup(node)[0]

		ages   = numpy.asarray([k.attrs.time for k in self.data.root.Z02])
		
		data = {}
		if _age in ages:
			#no interpolation, isochrone already in the file
			t = self.data.getNode('/Z'+str(metal)[2:]+'/a'+str(_age)) 
			for kn in t.colnames:
				data [kn] = t.col(kn) 
		else:
			#interpolate between isochrones
			d      = (ages - float(_age))**2
			a1, a2 = ages[numpy.argsort(d)[:2]]
			#print "Warning: Interpolation between %d and %d Myr" % (a1,a2)
			r = numpy.log10(_age/a1)/numpy.log10(a2/a1)

			t1 = self.data.getNode('/Z'+str(metal)[2:]+'/a'+str(int(a1)))	
			t2 = self.data.getNode('/Z'+str(metal)[2:]+'/a'+str(int(a2)))	

			stop = min(t1.nrows,t2.nrows)
		
			for kn in t1.colnames:
				y2 = t2.col(kn)[:stop]
				y1 = t1.col(kn)[:stop]
				data[kn] = y2 * r + y1 * (1.-r)
				del y1, y2

		#mass selection
		if masses != None:
			#masses are expected in logM for interpolation
			if masses.max() > 2.3:
				_m = numpy.log10(masses)
			else:
				_m = masses
			data_logM = data['logM'][:]
			for kn in data:
				data[kn] = interp(_m, data_logM, data[kn])
				
		table = mytables.Table(data, name='Isochrone from %s' % self.name)
		table.header['metal'] = metal
		table.header['time'] = _age * 1e6
		return table	
