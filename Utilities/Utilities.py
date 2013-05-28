import numpy

# =====================================
# Make Contour
# =====================================
def contour(x, y, F, base_directory='.', show_plot=False):
	return F[0, 0]*x**2 + F[1, 1]*y**2 + F[0, 1]*x*y + F[1, 0]*x*y

# ==============================
# define smoothing
# ==============================
def k_smooth_py(loglam,flux,vdisp):

	nlambda		  = len(loglam)
	pixsize		  = numpy.abs(numpy.log(10.)* 2.99792e5* (loglam[nlambda-1]-loglam[0])/float(nlambda))
	smoothing	  = float(vdisp)/float(pixsize)								 # pixels
	npix		  = long(4.0* numpy.ceil(smoothing))* long(2)+3
	klam		  = numpy.arange(float(npix))- float(npix- 1.)/ 2.
	kernel		  = numpy.exp(-0.5* (klam/ smoothing)**2)/numpy.sqrt(2.*numpy.pi)/smoothing
	kernel		  = kernel/sum(kernel)
	smoothed_spec = numpy.convolve(flux,kernel,mode='same')#,/edge_truncate)	python code can't do this yet.# "valid" mode gives too short

	return smoothed_spec


# ==============================
# Make Histogram
# ==============================
def make_histogram(array, nbins, mymin=-1, mymax =-1):

	# set z min/max
	if mymin == -1 : mymin  = numpy.min(array)
	if mymax == -1 : mymax  = numpy.max(array)

	# measure histogram
	# --- using bins=nbins (where nbins is a scalar number) makes the bins necessarily even
	hist,bin_edges = numpy.histogram(array,bins=nbins)

	# calculate bin centers from edges
	# --- there will be one less center value than the edges
	center = (bin_edges[:-1]+bin_edges[1:])/2

	# calculate bin widths, which are even
	width  = bin_edges[1]-bin_edges[0]

	return width, center, bin_edges, hist



# =========================================================================
# TOOLS for importing CSV files
# =========================================================================

def importArrayFromCSVFile(filename):
	# this imports the VISTA catalog (cmc_vista_des_...) and returns the contents as a NumPy data array

	#http://docs.python.org/2/library/csv.html#csv.Sniffer
	#dataArray = numpy.genfromtxt(filename, dtype='float32', comments='#')

	stringData = stringDataFromCSVFile(filename)

	# second level: prune comment lines
	stringDataPruned = pruneStringData(stringData)

	# third level: convert strings to floats
	# and collect them in a list of lists
	floatData = floatFromStringData(stringDataPruned)

	# fourth level: convert to NumPy array
	dataArray = numpy.array(floatData)

	#print filename, numpy.max((dataArray1 - dataArray + numpy.ones(dataArray.shape) / 1000000000)/dataArray)
	return dataArray

# needed subroutines

def stringDataFromCSVFile(filename):
	# create a CSV reader object
	csvreader = csv.reader(open(filename,'rb'),delimiter=' ',quotechar='|')

	# first level of extraction: strings
	stringData = []
	for row in csvreader:
		stringData.append(row)

	return stringData


def pruneStringData(stringData):
	# second level: discard comment lines and empty strings from whitespace
	stringDataPruned = []
	for row in stringData:
		if not (row[0].startswith('#')):
			prunedRow = []
			for entry in row:
				if (entry != ''):
					prunedRow.append(entry)
			stringDataPruned.append(prunedRow)

	return stringDataPruned




def floatFromStringData(stringDataPruned):
	#	floatData = []
	#	for row in stringDataPruned:
	#		stringRow = row[0]
	#		floatRow = [convertToFloat(s) for s in stringRow.split()]
	#		floatData.append(floatRow)
	floatData = []
	for row in stringDataPruned:
		floatRow = [convertToFloat(s) for s in row]
		floatData.append(floatRow)

	return floatData


def convertToFloat(string):
	if (string.startswith('*')):
		return 0
	else:
		return float(string)

# =========================================================================
# END OF CSV import tools
# =========================================================================
