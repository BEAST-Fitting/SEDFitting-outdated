import numpy

class Observations(object):

    def __init__(self, inputFile, distanceModulus=0., desc=None):
        """ Generate a data interface object """
        self.inputFile = inputFile
        self.filters   = None
        self.desc      = desc
        self.setDistanceModulus(distanceModulus)
        self.readData()
	self.badvalue  = None

    @property
    def nObs(self):
        return self.data.nrows

    def __len__(self):
	    return self.nObs

    def __call__(self):
        """ Calling the object will show info """
        print "Data read from %s " % self.inputFile
	if self.desc != None:
		print "Description: %s" % self.desc
        print "Number of records: %d" % self.nObs
        print ""
        print "Dataset contains:"

        for k in self.data.keys(): 
            print "\t %s" % k
	
	if self.filters == None:
		print "No filters set yet!"
	else:
		print "Using filters:", self.filters
    
    def __getitem__(self, *args, **kwargs):
        """ get item will generate a subsample """
	return self.data.__getitem__(*args, **kwargs)

    def keys(self):
        """ Returns dataset content names """
        return self.data.keys()

    def setDescription(self, txt):
	self.desc = txt

    def setDistanceModulus(self, val):
        """ Set the distance modulus to consider the dataset """
        self.distanceModulus = val
        self.distance = 10**( (val - 25.)/5. )
    
    def setDistance(self, val):
        """ Set observed object distance to X Megaparsecs 
            this will update also the distance Modulus
        """
        self.distance = val
        self.distanceModulus = 5. * log10( val*1e5 )
	
    def setBadValue(self, val):
	    self.badvalue = val

    def getFilters(self):
	return self.filters

    def setFilters(self, filters):
	self.filters = filters

    def getMags(self, num, filters):
	return numpy.array([ self.data[tt][num] - self.distanceModulus for tt in
		filters])

    def getErrors(self, num, filters):
	return numpy.array([ self.data[tt+'err'][num] for tt in filters])

    def getMask(self, num, filters, mags=None):
	if mags is None:
		mags = self.getMags(num, self.filters)

	if not self.badvalue is None:
		mask = ((mags >= self.badvalue) | numpy.isnan(mags) | numpy.isinf(mags) )
	else:
		mask = numpy.zeros(len(mags), dtype=bool)
	return mask

    def getObs(self, num=0):
        """ returns the dictionnary used during the analysis """
        assert ( not self.filters is None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = self.getErrors(num, self.filters)
	mask = self.getMask(num, self.filters)
        return mags, errs, mask

    def readData(self):
        """ read the dataset from the original source file """
	from tools import mytables
	self.data = mytables.load(self.inputFile)

    def __repr__(self):
	txt  = 'Observations: %s' % self.desc
	txt += '\n    distance = %g Mpc' % self.distance
	txt += '\n       |M-m| = %g mag' % self.distanceModulus
	txt += '\n    #objects = %d' % self.nObs
	txt += '\n%s' % object.__repr__(self)
	return txt





