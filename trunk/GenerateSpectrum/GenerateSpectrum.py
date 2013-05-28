# PURPOSE:
#   Reconstruct galaxy rest-frame spectrum given a fit
# CALLING SEQUENCE:
# INPUTS:
#   coeffs - [NT, NGALS] coefficients
# OUTPUTS:
#   loglam - [NL, NGALS] wavelengths
#   flux - [NL, NGALS] fluxes (ergs cm-2 s-1 A-1)
# OPTIONAL INPUTS:
# OPTIONAL OUTPUTS:
#   nt - total number of templates
#   mass, metallicity - properties of template fit;
#					   mass is current stellar mass and is in units of
#					   1 solar mass / (D/10pc)**2
#   b300 - star-formation within last 300Myrs relative to average
#		  star-formation rate
#   b1000 - star-formation within last 1Gyrs relative to average
#		   star-formation rate
#   lspec = spectrum of lines [0 , 0 , 0 , SOMETHING ,  0 0 0 0 ]
#
# implementation of stellar pieces
# 	optional outputs typically for galaxy evolution, but not cosmo
#	mass		= sum(tmremain* coeffs)
#	b300		= sum(tmass300* coeffs)/ sum(tmass* coeffs)  #stars gen in last 300myr
#	b1000		= sum(tmass1000* coeffs)/ sum(tmass* coeffs)
#	metallicity	= sum(tmremain* tmetallicity* coeffs)/ mass
#
# COMMENTS:
#   If coeffs are standard, returns units of erg/cm**2/s/A
#   If vdisp, /nolines, and /noextinct
#   If lines are included, they are always smoothed at 300 km/s vdisp
#   Bases fit on file:
#	  $KCORRECT_DIR/data/templates/k_nmf_derived.[vname].fits

# Questions to be answered by Cunha, Busha
#  do we smooth the continuum, or use the default of 300 km/s?  Will usually want vdisp = 100 for DESpec
#  Note that this assumes that smoothing is velocity dependant, not wavelength dependant.  DESpec could be constant in wavelength, so we might need to add this capability
#  also, do we include lines or extinction?
#------------------------------------------------------------------------------
import sys
sys.path.append('../../Wrapper/')
import archive
sys.path.append('../../Plotters/')
import Plotters
sys.path.append('../../Utilities/')
import Utilities

import numpy
from scipy.signal import convolve


# global variable for data bank file
databank_file = '../../../data/data_bank.h5'

class GenerateSpectrum:


	def __init__(self):
		print
		print
		print "========================================"
		print "   Reconstruct Pure Galaxy Spectrum  "
		print "========================================"

		#=======================================================================================
		# Input
		#=======================================================================================
		print "... Import data and params"
		with archive.archive(databank_file, 'r') as ar:

			# observation choices
			self.nolines	 = ar['/GenerateSpectrum/nolines']  	# user inputs
			self.noextinct	 = ar['/GenerateSpectrum/noextinct']
			self.vdisp		 = ar['/GenerateSpectrum/vdisp']

			# templates
			self.loglam 		  = numpy.array(ar['/spectral_templates/log10_wavelength'])
			self.tspec_v0 		  = numpy.array(ar['/spectral_templates/template_spectrum_v0']) # not smoothed
			self.tspec_v0_nl 	  = numpy.array(ar['/spectral_templates/template_spectrum_v0_nolines'])
			self.tspec_v0_nd 	  = numpy.array(ar['/spectral_templates/template_spectrum_v0_nodust'])
			self.tspec_v0_nd_nl	  = numpy.array(ar['/spectral_templates/template_spectrum_v0_nodust_nolines'])
			self.tspec_v300 	  = numpy.array(ar['/spectral_templates/template_spectrum_v300']) # smoothed to 300kms(res)
			self.tspec_v300_nl 	  = numpy.array(ar['/spectral_templates/template_spectrum_v300_nolines'])
			self.tspec_v300_nd 	  = numpy.array(ar['/spectral_templates/template_spectrum_v300_nodust'])
			self.tspec_v300_nd_nl = numpy.array(ar['/spectral_templates/template_spectrum_v300_nodust_nolines'])
			self.lspec_v300 	  = numpy.array(ar['/spectral_templates/lspec_v300'])

			# galaxy information
			self.coeffs		 = ar['/gal/coeffs']
			self.z_true		 = ar['/gal/z_true']
			self.fiber_selection_flag = ar['/gal/fiber_selection_flag']

		print "... finished importing data"


		# =====================================================
		# QA Test -- Incoming
		Ngal  = len(self.fiber_selection_flag)	 # count galaxies
		Nwave = len(self.tspec_v0[0,:])	 # count wavelengths

		fiber_selection_flag_shape_expected = (Ngal,)
		check = (self.fiber_selection_flag.shape == fiber_selection_flag_shape_expected)
		if check :
			print "*** Incoming QA checksum PASSED ***"
		else:
			print "XX CHECKSUM FAILED XX"
			print "fiber selection flag shape not compliant withe expectation"
			sys.exit(1)
	# =====================================================



	def clean_spectrum(self, galaxy_index):

		# set default velocity dispersion
		sqrt_vdisp		 = numpy.sqrt(self.vdisp)
		coeffs_temporary = self.coeffs[galaxy_index]
		z_temporary		 = self.z_true[galaxy_index]


		# Create Spectra from options
		if self.nolines == 1 and self.noextinct == 1 :
			flux_temporary = self.tspec_v0_nl_nd* coeffs_temporary
			flux_temporary = Utilities.k_smooth_py(self.loglam,flux_temporary,sqrt_vdisp)

		if self.nolines == 1 and self.noextinct == 0 :
			flux_temporary = self.tspec_v0_nl* coeffs_temporary
			flux_temporary = Utilities.k_smooth_py(self.loglam,flux_temporary,sqrt_vdisp)

		if self.nolines == 0 and self.noextinct == 1 :
			flux_temporary  = self.tspec_v0_nd* coeffs_temporary
			flux_temporary  = Utilities.k_smooth_py(self.loglam,flux_temporary,sqrt_vdisp)
			flux_temporary += self.lspec_v300* coeffs_temporary


		#  should be the default option with vdisp = 100
		if self.nolines == 0 and self.noextinct == 0 :

			# redshift evolve the coefficients
			coeffs_temporary_z = coeffs_temporary* (1.+ z_temporary)**1.8

			# set up the flux and make the first spectrum combination
			flux_temporary = numpy.zeros(self.tspec_v0[0,:].shape)
			for j in range(coeffs_temporary.shape[0]) :
				flux_temporary += coeffs_temporary[j]*self.tspec_v0[j,:]

			# smooth spectrum for the given velocity dispersion
			flux_temporary	 = numpy.array(Utilities.k_smooth_py(self.loglam,flux_temporary,sqrt_vdisp))

			for j in range(len(coeffs_temporary)):
				flux_temporary += coeffs_temporary_z[j]*self.lspec_v300[j,:]

		# flux dimming and redshift stretching via redshift
		flux_temporary = flux_temporary/ (1. + z_temporary)


		#=======================================================================================
		# Redshift galaxy sed's.
		#	   - rest-frame flux assumes galaxy is 10pc away.
		#	   - doesn't matter when using coeffs calculated from
		#	   - kcorrect, since the code already applies the 1/r^2 dimming.
		#=======================================================================================
		loglam_temp = self.loglam+ numpy.log10(1. + z_temporary)

		# re-scale lambda (wavelengths)
		lam_new = 10**loglam_temp

		return flux_temporary, lam_new


	#=======================================================================================
	# Loop over galaxies create all spectra
	#=======================================================================================
	def loop_over_gal(self):

		#=======================================================================================
		#  LOOP over galaxies	 # default add both ... 0 for both # 1: don't do it
		#=======================================================================================
		check_fiber_sel = numpy.where(self.fiber_selection_flag == True)
		self.coeffs = self.coeffs[check_fiber_sel,:]
		self.coeffs = numpy.reshape(self.coeffs,self.coeffs.shape[1:])
		self.z_true = self.z_true[check_fiber_sel]

		flux_new		= []
		lam_new			= []
		galax_index_new = []

		Ngal_new = len(check_fiber_sel[0])
		print "... Loop over galaxies to create spectral templatesi [this number has been cut]", Ngal_new
		for i in range(Ngal_new):

			# calculate spectrum from coefficients
			flux_temporary, lam_temp  = self.clean_spectrum(i)

			# Add temporary to output/new
			flux_new.append(numpy.array(flux_temporary))
			lam_new.append(numpy.array(lam_temp))


		#---------------------------------------------------------------------------------------
		# END LOOP over galaxies
		#---------------------------------------------------------------------------------------

		# make output into numpy arrays
		flux_new = numpy.array(flux_new)
		lam_new  = numpy.array(lam_new)

		# DIAGNOSTIC
		print "... plotting diagnostics"
		Plotters.plot_spectrum(lam_new[0,:], flux_new[0,:], xmin=5000, xmax=10000,fig_number=0, plot_number=0, title='Galaxy Spectrum', base_directory='.', filename_addendum="", show_plot=show_plot)


		print "... Loop finished, ngal with spectra = ", len(flux_new)
		print "... shape of output arrays:"
		print "... ... wavelength",lam_new.shape
		print "... ... flux		 ",flux_new.shape


		# Output
		print "... output to databank"
		with archive.archive(databank_file, 'a') as ar:
			ar['/gal/wavelength'] = lam_new
			ar['/gal/flux']		  = flux_new




#=======================================================================================
# Main
#=======================================================================================
def main():

	s = GenerateSpectrum()
	s.loop_over_gal()

	sys.exit(0)



#=======================================================================================
if __name__ == '__main__':
	main()
