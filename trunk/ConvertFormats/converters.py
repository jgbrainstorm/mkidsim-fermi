# 17.1.2013
# Ben Hambrecht
# Cosmology Research Group, ETH
# hambrecht@phys.ethz.ch

import sys
sys.path.append('../Wrapper/')
import archive
import ConfigParser
import time
import numpy
import pyfits
import csv
import asciitable
import os
import shutil
import string
import re
sys.path.append('../Wrapper/')
import archive
from archive import nodeList




# ==================================================
#
# ==================================================
def bcc_hdf_converter(observed_cat,truth_cat,outputfile,max_gals, min_ra=[], min_dec=[]):

	# OBSERVED CAT #

	# open FITS file
	catalog = pyfits.open(observed_cat)

	# extracting data table
	catalogData = catalog[1].data

	# print catalog[1].columns.names

	#extracting columns (http://www.slac.stanford.edu/~risa/des_bcc/bcc_v0.5_tags.html)
	galaxy_index = catalogData.field('id') # 'index'?
	magnitude_U  = catalogData.field('mag_u')
	magnitude_G  = catalogData.field('mag_g')
	magnitude_R  = catalogData.field('mag_r')
	magnitude_I  = catalogData.field('mag_i')
	magnitude_Z  = catalogData.field('mag_z')
	magnitude_Y  = catalogData.field('mag_y')
	redshift_photo_gaussian = catalogData.field('photoz_gaussian')
	ra_true  = catalogData.field('ra')
	dec_true = catalogData.field('dec')

	# ---------------------------
	# set up alternative mininum for RA/DEC in order to shift galaxy coordinates to the region of sky we look at (J. Forero-Romero)
	if(min_ra and min_dec):

		# obtain minima of original ra/dec in truth tables
		min_ra_true = numpy.amin(ra_true)
		min_dec_true = numpy.amin(dec_true)
		deg_to_radian = numpy.pi/180.0

		#recompute RA from the expected DEC positions
		ra_true = ra_true * (numpy.cos(deg_to_radian * dec_true))/(numpy.cos(deg_to_radian * (dec_true - min_dec_true + min_dec)))

		#make the final shift in DEC
		dec_true = dec_true - min_dec_true + min_dec

		#make the final shift in RA
		ra_true = (ra_true - numpy.amin(ra_true) + min_ra)


	print 'Final min max ra', numpy.amin(ra_true), numpy.amax(ra_true)
	print 'Final min max dec', numpy.amin(dec_true), numpy.amax(dec_true)
	# ---------------------------

	# close FITS file
	catalog.close()

	number_of_gals = len(galaxy_index)
	numpy.random.seed(seed=99999)

	if (max_gals == 0):
		max_gals = number_of_gals
		random_selection_flag = (numpy.random.random_sample((number_of_gals,)) > -1)
	else:
		random_selection_flag = (numpy.random.random_sample((number_of_gals,))*number_of_gals < max_gals)


	# write to HDF5 file
	with archive.archive(outputfile,'w') as ar:
		ar['/gal/galaxy_index'] = galaxy_index[random_selection_flag]
		ar['/gal/magnitude_u'] = magnitude_U[random_selection_flag]
		ar['/gal/magnitude_g'] = magnitude_G[random_selection_flag]
		ar['/gal/magnitude_r'] = magnitude_R[random_selection_flag]
		ar['/gal/magnitude_i'] = magnitude_I[random_selection_flag]
		ar['/gal/magnitude_z'] = magnitude_Z[random_selection_flag]
		ar['/gal/magnitude_y'] = magnitude_Y[random_selection_flag]
		ar['/gal/redshift_photometric'] = redshift_photo_gaussian[random_selection_flag]
		ar['/gal/ra_true']  = ra_true[random_selection_flag]
		ar['/gal/dec_true'] = dec_true[random_selection_flag]



#		ar['/gal/min_ra_true'] =
#		ar['/gal/max_ra_true'] =
#		ar['/gal/min_dec_true'] =
#		ar['/gal/max_dec_true'] =

	# TRUTH CAT #

	# open FITS file
	catalog = pyfits.open(truth_cat)

	# extracting data table
	catalogData = catalog[1].data

	#extracting columns (http://www.slac.stanford.edu/~risa/des_bcc/bcc_v0.5_tags.html)
	coeffs = numpy.array(catalogData.field('coeffs'))
	z_true = numpy.array(catalogData.field('z'))

	# close FITS file
	catalog.close()

	# write to HDF5 file
	with archive.archive(outputfile,'a') as ar:
		ar['/gal/coeffs'] = coeffs[random_selection_flag]
		ar['/gal/z_true'] = z_true[random_selection_flag]


# ==================================================
#
# ==================================================
def entriesInBCCCatalog(catalog):
	# What kind of quantities does the catalog contain for each galaxy?
	return catalog[1].columns.names





# ==================================================
#
# ==================================================
def cosmos_hdf_converter(inputfile,outputfile):

	dataArray = importArrayFromCSVFile(inputfile)

	# careful: index here is shifted by 1 wrt the index in the COSMOS header
	galaxy_index = dataArray[:,0]

	photometric_redshift = dataArray[:,1]
	spectroscopic_redshift = dataArray[:,31]

	magnitude_g = dataArray[:,14]
	magnitude_r = dataArray[:,15]
	magnitude_i = dataArray[:,16]
	magnitude_z = dataArray[:,17]
	magnitude_Y = dataArray[:,18]
	magnitude_J = dataArray[:,19]
	magnitude_H = dataArray[:,20]
	magnitude_K = dataArray[:,21]
	magnitude_g_error = dataArray[:,22]
	magnitude_r_error = dataArray[:,23]
	magnitude_i_error = dataArray[:,24]
	magnitude_z_error = dataArray[:,25]
	magnitude_Y_error = dataArray[:,26]
	magnitude_J_error = dataArray[:,27]
	magnitude_H_error = dataArray[:,28]
	magnitude_K_error = dataArray[:,29]

	#???ra = dataArray[:,XXX]  # Stephanie gave us a new file we just need to import
	#???dec = dataArray[:,XXX]

	# write to HDF5 file
	with archive.archive(outputfile,'a') as ar:
		ar['/gal/galaxy_index'] = galaxy_index

		ar['/gal/redshift_photometric'] = photometric_redshift
		ar['/gal/redshift_spectroscopic'] = spectroscopic_redshift

		ar['/gal/magnitude_g'] = magnitude_g
		ar['/gal/magnitude_r'] = magnitude_r
		ar['/gal/magnitude_i'] = magnitude_i
		ar['/gal/magnitude_z'] = magnitude_z
		ar['/gal/magnitude_y'] = magnitude_Y
		ar['/gal/magnitude_j'] = magnitude_J
		ar['/gal/magnitude_h'] = magnitude_H
		ar['/gal/magnitude_k'] = magnitude_K

		ar['/gal/magnitude_error_g'] = magnitude_g_error
		ar['/gal/magnitude_error_r'] = magnitude_r_error
		ar['/gal/magnitude_error_i'] = magnitude_i_error
		ar['/gal/magnitude_error_z'] = magnitude_z_error
		ar['/gal/magnitude_error_y'] = magnitude_Y_error
		ar['/gal/magnitude_error_j'] = magnitude_J_error
		ar['/gal/magnitude_error_h'] = magnitude_H_error
		ar['/gal/magnitude_error_k'] = magnitude_K_error

# =========================================================================
def cosmos_insert_ra_dec(inputfile,outputfile):
	dataArray = importArrayFromCSVFile(inputfile)

	right_ascension = dataArray[:,1]
	declination = dataArray[:,2]

	# append to HDF file
	with archive.archive(outputfile,'a') as ar:
		ar['/gal/ra_true']  = right_ascension
		ar['/gal/dec_true'] = declination





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





def importEmissionLineSensitivityFromCSVFile(inputfile,outputfile):

	dataArray = importArrayFromCSVFile(inputfile)

	wavelength = dataArray[:,0] * 1e-10 # conversion Angstroem to meter
	signal_to_noise_at_3sigma = dataArray[:,1]
	signal_to_noise_at_5sigma = dataArray[:,2]

	# write to HDF5 file
	with archive.archive(outputfile,'a') as ar:
		ar['/instrument/sensitivity/emission_line_wavelengths'] = wavelength
		#		ar['/gal/wavelength/unit'] = 'meter'
		ar['/instrument/sensitivity/signal_to_noise_sigma3'] = signal_to_noise_at_3sigma
		ar['/instrument/sensitivity/signal_to_noise_sigma5'] = signal_to_noise_at_5sigma
		ar['/instrument/sensitivity/exposure_time'] = 3000



def importContinuumSensitivityFromCSVFile(inputfile,outputfile):

	dataArray = importArrayFromCSVFile(inputfile)

	wavelength = dataArray[:,0] * 1e-10 # conversion Angstroem to meter
	signal_to_noise_continuum = dataArray[:,1]
	signal_to_noise_emission_line = dataArray[:,2]

	# write to HDF5 file
	with archive.archive(outputfile,'a') as ar:
		ar['/instrument/sensitivity/wavelength_grid'] 			= wavelength
		ar['/instrument/sensitivity/signal_to_noise_continuum']	 = signal_to_noise_continuum
		ar['/instrument/sensitivity/signal_to_noise_emission_line'] = signal_to_noise_emission_line
#		ar['/instrument/sensitivity/exposure_time'] = 3000




def importTrainingCatalog(inputfile,outputfile):

	dataArray = importArrayFromCSVFile(inputfile)

	redshift_photometric = dataArray[:,1]
	galaxy_type = dataArray[:,2]

	magnitude_g = dataArray[:,5]
	magnitude_r = dataArray[:,8]
	magnitude_i = dataArray[:,11]
	magnitude_z = dataArray[:,14]
	magnitude_y = dataArray[:,17]
	magnitude_j = dataArray[:,20]
	magnitude_h = dataArray[:,23]
	magnitude_k = dataArray[:,26]

	#	wavelength_Ly = dataArray[:,29] * 1e-10 # conversion Angstroem to meter
	#	wavelength_OII = dataArray[:,32] * 1e-10
	#	wavelength_Hb = dataArray[:,35] * 1e-10
	#	wavelength_OIIIa = dataArray[:,38] * 1e-10
	#	wavelength_OIIIb = dataArray[:,41] * 1e-10
	#	wavelength_Ha = dataArray[:,44] * 1e-10



	# write to HDF5 file
	with archive.archive(outputfile,'a') as ar:
		ar['/target_selection_training/gal/redshift_photometric'] = redshift_photometric
		ar['/target_selection_training/gal/galaxy_type'] = galaxy_type
		ar['/target_selection_training/gal/magnitude_g'] = magnitude_g
		ar['/target_selection_training/gal/magnitude_r'] = magnitude_r
		ar['/target_selection_training/gal/magnitude_i'] = magnitude_i
		ar['/target_selection_training/gal/magnitude_z'] = magnitude_z
		ar['/target_selection_training/gal/magnitude_y'] = magnitude_y
		ar['/target_selection_training/gal/magnitude_j'] = magnitude_j
		ar['/target_selection_training/gal/magnitude_h'] = magnitude_h
		ar['/target_selection_training/gal/magnitude_k'] = magnitude_k
#		ar['/target_selection_training/gal/emission_lines/wavelength_Ly'] = wavelength_Ly
#		ar['/target_selection_training/gal/emission_lines/wavelength_OII'] = wavelength_OII
#		ar['/target_selection_training/gal/emission_lines/wavelength_Hb'] = wavelength_Hb
#		ar['/target_selection_training/gal/emission_lines/wavelength_OIIIa'] = wavelength_OIIIa
#		ar['/target_selection_training/gal/emission_lines/wavelength_OIIIb'] = wavelength_OIIIb
#		ar['/target_selection_training/gal/emission_lines/wavelength_Ha'] = wavelength_Ha




# to be used for converting spectral templates for Sim_Gal_clean
def import_convert_spectra_templates(inputfile,outputfile):

	#input files
	#inputfile	=  'k_nmf_derived.default.fits' 	# make this file name a parameter in one of the ini files
	#inputfile	=  'k_nmf_derived.newdefault.fits'

	# open file
	hdunum = 0
	hdu = pyfits.open(inputfile)[hdunum]

	# header
	## read in header
	header = hdu.header

	## read in header element
	nt	= header['nt']
	#print '... number of templates', nt

	# read in table elements
	loglam	 		= numpy.array(numpy.log10(pyfits.open(inputfile)[11].data))
	tspec_v0		= numpy.array(pyfits.open(inputfile)[1].data)
	tspec_v0_nl		= numpy.array(pyfits.open(inputfile)[2].data)
	tspec_v0_nd		= numpy.array(pyfits.open(inputfile)[3].data)
	tspec_v0_nd_nl	= numpy.array(pyfits.open(inputfile)[4].data)
	tspec_v300		= numpy.array(pyfits.open(inputfile)[5].data)
	tspec_v300_nl	= numpy.array(pyfits.open(inputfile)[6].data)
	tspec_v300_nd	= numpy.array(pyfits.open(inputfile)[7].data)
	tspec_v300_nd_nl= numpy.array(pyfits.open(inputfile)[8].data)
	lspec_v300		= numpy.array(pyfits.open(inputfile)[9].data)
	#	tmass			= pyfits.open(inputfile)[17].data
	#	tmetallicity		= pyfits.open(inputfile)[18].data
	#	tmass300		= pyfits.open(inputfile)[19].data
	#	tmass1000		= pyfits.open(inputfile)[20].data
	#	tmremain		= pyfits.open(inputfile)[24].data


	# write to databank
	with archive.archive(outputfile,'a') as ar:
		ar['/spectral_templates/num_spectral_templates'] 				= nt		# number of spectral spectral_templates; nt
		ar['/spectral_templates/log10_wavelength'] 						= loglam
		ar['/spectral_templates/template_spectrum_v0'] 					= tspec_v0			  # not smoothed
		ar['/spectral_templates/template_spectrum_v0_nolines'] 			= tspec_v0_nl
		ar['/spectral_templates/template_spectrum_v0_nodust'] 			= tspec_v0_nd
		ar['/spectral_templates/template_spectrum_v0_nodust_nolines'] 	= tspec_v0_nd_nl
		ar['/spectral_templates/template_spectrum_v300']	 			= tspec_v300			# smoothed to 300kms(res)
		ar['/spectral_templates/template_spectrum_v300_nolines'] 		= tspec_v300_nl
		ar['/spectral_templates/template_spectrum_v300_nodust'] 		= tspec_v300_nd
		ar['/spectral_templates/template_spectrum_v300_nodust_nolines'] = tspec_v300_nd_nl
		ar['/spectral_templates/lspec_v300'] 							= lspec_v300
#		ar['/spectral_templates/template_mass'] =  tmass				 # stellar stuff
#		ar['/spectral_templates/template_metallicity'] =  tmetallicity
#		ar['/spectral_templates/template_mass300'] =  tmass300
#		ar['/spectral_templates/template_mass1000'] =  tmass1000
#		ar['/spectral_templates/template_mass_remain'] =  tmremain


# read in noise files for making noisy observation in Sim_Gal_Spec_observed
def convert_sky_bkg(inputfile,outputfile): # from skybg_50_10

	table = asciitable.read(inputfile)
	wave = []
	power = []
	for i in range(len(table)):
		wave.append(table[i][0])
		power.append(table[i][1])

	with archive.archive(outputfile,'a') as ar:
		ar['/sky/background/wave'] = numpy.array(wave)
		ar['/sky/background/power'] = numpy.array(power)

# read in noise files for making noisy observation in Sim_Gal_Spec_observed
def convert_extinction(inputfile,outputfile): # from palomar

	table = asciitable.read(inputfile)
	abswave = []
	abstemp = []
	for i in range(len(table)):
		abswave.append(table[i][0])
		abstemp.append(table[i][1])

	with archive.archive(outputfile,'a') as ar:
		ar['/sky/extinction/wave'] = numpy.array(abswave)
		ar['/sky/extinction/power'] = numpy.array(abstemp)

# =====================================================
def main():

	data_bank_location = '../../data/data_bank.h5'

	with archive.archive(data_bank_location,'r') as ar:
		nb_galaxies_sample = ar['/ConvertFormats/nb_galaxies_imported']
		truth_catalog_location = ar['/ConvertFormats/truth_catalog_filename']
		observed_catalog_location = ar['/ConvertFormats/observed_catalog_filename']
		spectral_templates_location = ar['/ConvertFormats/spectral_templates_filename']
		sky_background_location = ar['/ConvertFormats/sky_background_filename']
		extinction_location = ar['/ConvertFormats/extinction_filename']
		throughput_parameter_location = ar['/ConvertFormats/throughput_parameters_filename']
		#recreate_training = ar['/Convert/recreate_training_only']
	recreate_training = True


	print
	print
	print "============================="
	print "		 Convert				"
	print "============================="



	with archive.archive(data_bank_location,'r') as ar:
		min_ra = ar['/ConvertFormats/min_right_ascension']
		min_dec = ar['/ConvertFormats/min_declination']
		import_training_catalog = ar['/ConvertFormats/import_training_catalog']

	# remove old copies
	if (os.path.isfile(data_bank_location)):
		os.remove(data_bank_location)


	#print "...importing BCC"
	# convert BCC catalog (FITS) to HDF5
	# the keywords for accessing the columns in the catalog change, see converters.py for details!
	bcc_hdf_converter(observed_cat=observed_catalog_location,truth_cat=truth_catalog_location,outputfile=data_bank_location,max_gals=nb_galaxies_sample, min_ra=min_ra, min_dec=min_dec) # this is the first so it has 'w', all others have 'a'
	print "... ... done [importing BCC]"

	# read out emission lines from WISE catalog
	if (import_training_catalog):
		#print "...importing Training Cat"
		importTrainingCatalog('../../data/COSMOS/CMC_VISTA_DES_jun11.cat',data_bank_location)
		print "... ... done [importing Training Cat]"

	#print "... importing spectral templates"
	import_convert_spectra_templates(spectral_templates_location,outputfile=data_bank_location)
	print "... ... done [importing spectral templates]"

	# read in noise files for making noisy observation in Sim_Gal_Spec_observed
	#print "... importing Sky"
	convert_sky_bkg(inputfile=sky_background_location,outputfile=data_bank_location)

	# read in noise files for making noisy observation in Sim_Gal_Spec_observed
	convert_extinction(inputfile=extinction_location,outputfile=data_bank_location)
	print "... ... done [importing Sky]"

	#import telescope information
	#print "... importing telescope throughput"
	import_throughput(throughput_parameter_location,outputfile=data_bank_location)
	print "... ... done [importing telescope throughput]"






# =====================================================
if __name__ == '__main__':
	main()













