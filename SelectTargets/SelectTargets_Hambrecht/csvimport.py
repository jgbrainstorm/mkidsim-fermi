import csv
import numpy



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
    #    floatData = []
    #    for row in stringDataPruned:
    #        stringRow = row[0]
    #        floatRow = [convertToFloat(s) for s in stringRow.split()]
    #        floatData.append(floatRow)
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
