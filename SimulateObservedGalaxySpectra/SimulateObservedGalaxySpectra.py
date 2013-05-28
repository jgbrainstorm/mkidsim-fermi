import sys
sys.path.append('../../Wrapper/')
sys.path.append('../../Plotters/')
sys.path.append('kcorrect_python_rewrite/')
sys.path.append('../../Utilities/')
import Utilities

import archive
import Plotters
import k_binspec_py as k_binspec

import subprocess
import numpy
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import random
from math import exp, sqrt
#from kcorrect import k_binspec
import pylab
import time

# set databank file
data_bank = '../../../data/data_bank.h5'
show_plot=False


# tell MSR to use method Cunha
with archive.archive(data_bank,'a') as ar:
	ar['/MeasureSpectroscopicRedshift/method'] = 'Cunha'



class SimulateObservedGalaxySpectra:
	
	def __init__(self):
		
		#=======================================================================================
		# Input
		#=======================================================================================
		print
		print
		print '========================='
		print 'Entering Noise Generator'
		print '========================='
		
		with archive.archive(data_bank,'r') as ar:
			
			# instrument information
			self.exposure  = 20.*60.				#exposure				= ar['/instrument/sensitivity/exposure_time']  # exposure = exposure time in seconds
			self.readnoise = 5.					  #readnoise 				= ar['/instrument/sensitivity/read_noise']			# add to param files
			#			resnum = .25							  #resnum 				= ar['/instrument/sensitivity/resolution']
			
			self.Nwave     = ar['/SimulateObservedGalaxySpectra/nb_pixels']
			wavelength_range = ar['/SimulateObservedGalaxySpectra/wavelength_range']
			wavelength_min = wavelength_range[0]
			wavelength_max = wavelength_range[1]
			
			fiber_aperture_radius = numpy.sqrt(2.5/numpy.pi)	#fiber_aperture_radius 	= ar['/Fiber_Allocation/fiber_radius']	  # arcsec
			self.altitude = 2635.				 #altitude 				= ar['/instrument/sensitivity/altitude']			 # meters  ...add to param files
			self.worsening = 1					   #worsening 				= ar['/instrument/sensitivity/worsening']			# add to param files
			
			self.z					= ar['/gal/z_true']
			
			tile_id_gal		        = ar['/gal/tile_ID']
			airmass_tile			= ar['/tiling/airmass']
			tile_id_tile		    = ar['/tiling/tile_ID']
			self.air_mass			= ar['/gal/airmass']
			
			# general constants
			self.hp=6.62607e-27			  #ergs.secs
			self.c=2.9979e18				 #Ang/secs
			self.seed_constant = 99999
			
			# Throughput
			self.wavetrans   = numpy.array(ar['/Throughput/parameters/wavelength'])
			self.trans		= numpy.array(ar['/Throughput/throughput'])  #!!! what should be the units here?
			
			# sky data
			self.wave   		= ar['/sky/background/wave']
			self.atmtemp		= ar['/sky/background/power']
			self.abswave 		= ar['/sky/extinction/wave']
			self.abstemp 		= ar['/sky/extinction/power']
			
			# basic modifications to input data
			self.atmtemp	= numpy.array(self.atmtemp, dtype='float64')
			self.wave		= numpy.array(self.wave, dtype='float64')
			self.abstemp	= numpy.array(self.abstemp)
			self.abswave	= numpy.array(self.abswave)
			
			
			#=======================================================================================
			# obtain airmasses for galaxies
			#=======================================================================================
            #for i in range(Number_Tiles):
            #    match_tile = numpy.where(tile_id_gal == tile_id_tile[i])[0]
            #    self.airmass[match_tile] = airmass_tile[i]
			
			#=======================================================================================
			# Set constants
			#=======================================================================================
			# !!! to modify to access tile information to link to where galaxy is observed
			self.num_readout = 1
			#Nwave		= 557*resnum		# 557 = number of pixels
			#			self.Nwave		= 1000*resnum		# 557 = number of pixels
			resnum = float(self.Nwave)/1000
			dispersion  = 2.5		# def: dispersion = (pixsize or width ) / resolution; 2.5 is when most of the information comes # in a limit
			#			dellam		= 7.14/resnum		#R is not being used. Only dellam matters now. = resolution/dispersion
			dellam = float(wavelength_max - wavelength_min)/self.Nwave
			
			print "Resolution:", int(self.Nwave), "pixels"
			
			
			#=======================================================================================
			# Aperture
			#=======================================================================================
			self.aperture_area 	= numpy.pi*fiber_aperture_radius**2
			self.Area			= numpy.pi*4**2
			self.Areacm		= self.Area*10000.		  #converting from m^2 to cm^2
			
			
			#=======================================================================================
			# SKY BACKGROUND
			# 	- Convert sky photons to photons/Ang.
			#	- in skybg_50_10.dat the units are photons/s/nm/arcsec^2/m^2
			#=======================================================================================
			self.atm	 = self.atmtemp* self.worsening* self.exposure* self.aperture_area* self.Area/10. #because 10ang=1 nm
			self.wave   *= 10. # converting angstroms
			
			
			#====================================
			# smooth spectra to SURVEY resolution
			# I think I shouldn't smooth over atmospheric absorption since
			# smoothing only happens after the photons have been removed from atmosphere.
			#====================================
			self.atm 	= Utilities.k_smooth_py(numpy.log10(self.wave), self.atm, 50)
			self.atm	= numpy.array(self.atm, dtype='float64')
			
			
			#=======================================================================================
			# TRANSMISSION information
			#=======================================================================================
			self.wavetrans   *= 10.  # converting angstroms
			#trans	   *= 0.01
			
			
			#=======================================================================================
			# Interpolate wavelengths onto grid with resolution
			#	1) make grid
			#	2) interpolate counts/fluxes onto grid
			#	3) Note that delta_lam will affect number of photons counted.
			#		***remember to multiply by delta_lam.***
			#=======================================================================================
			self.lambda_survey	 = numpy.zeros(self.Nwave, dtype='float64')
			binwidth		 = numpy.zeros(self.Nwave)
			self.lambda_survey[0] = wavelength_min
			
			for i in numpy.arange(1,self.Nwave):
				self.lambda_survey[i] = self.lambda_survey[i-1]+dellam
			
			for i in numpy.arange(1,self.Nwave-1):
				self.lambda_survey[i] = self.lambda_survey[i] + (self.lambda_survey[i+1]-self.lambda_survey[i])/2.
				binwidth[i] 	 = self.lambda_survey[i+1] - self.lambda_survey[i]
			
			atmtemp_rebinned = numpy.zeros(self.lambda_survey.shape)
			k_binspec.k_binspec(self.wave,self.atmtemp,self.lambda_survey,atmtemp_rebinned)
			self.trans_interp = numpy.interp(self.lambda_survey,self.wavetrans,self.trans)
			
			
			self.A = self.hp*self.c/(self.lambda_survey*self.trans_interp)
			self.A[numpy.isnan(self.A)] = 0.
			
			
			
			self.B = self.worsening*self.exposure*self.aperture_area*self.Area/10*atmtemp_rebinned*self.trans_interp
			
			self.C = self.exposure*self.Areacm/(self.hp*self.c)*self.lambda_survey*self.trans_interp
			
			abstemp_interp = numpy.interp(self.lambda_survey,self.abswave,self.abstemp)
			self.D = abstemp_interp*numpy.exp((1700.-self.altitude)/7000.)
			
			self.E = self.num_readout*self.readnoise**2
	
	
	
	def observed_spectrum(self, lam, flux, galaxy_index):
		
		air_mass_temp = self.air_mass[galaxy_index]
		#			#============================================
		#			# Modify ATMOSPHERIC ABSORPTION by location and airmass of observation
		#			#============================================
		#			absatm = self.abstemp*air_mass_temp*exp((1700.-self.altitude)/7000.)
		#			absatm2 = absatm/air_mass_temp
		##			print numpy.min(absatm2), numpy.max(absatm2)
		#			where_bad = numpy.where(absatm > 1)
		#			absatm[numpy.where(absatm > 1)] = 1.			 # file is bad ,...just make sure file is not bad
		#
		#			Vtrans = self.trans_interp
		#		
		#			#print '================================='
		#			#print 'check nan, range'
		#			#=======================================================================================
		#			# Convert sed flux to photon counts.
		#			# 	- sed flux is in ergs/cm^2/s/Ang.
		#			# 	- photon counts is in photons/Ang.
		#			#=======================================================================================
		#			sedphot = self.exposure* flux* lam* self.Areacm/ (self.hp*self.c) 	# photons/Ang   #approximate - correct is to integrate
		#			sedflux = self.exposure* flux* self.Areacm				# ergs/Ang
		#
		#
		#			#=======================================================================================
		#			#	Bin spectra to survey's resolutions by integrating over pixels
		#			#		to convert flux density to flux.
		#			#=======================================================================================
		#
		#			# cut wavelengths to prevent trying to extrapolate during rebinning/interpolation
		#			cut_indexes = numpy.where(numpy.logical_and((lam > self.wave[0]),(lam < self.wave[-1])))
		#			lam = lam[cut_indexes]
		#			sedphot = sedphot[cut_indexes]
		#			sedflux = sedflux[cut_indexes]
		#
		#
		#			# create new arrays
		#			Vatm =  numpy.zeros(len(self.lambda_survey), dtype='float64')
		#			Vsed =  numpy.zeros(len(self.lambda_survey), dtype='float64')
		#			Vflux = numpy.zeros(len(self.lambda_survey), dtype='float64')
		flux_rebinned = numpy.zeros(len(self.lambda_survey), dtype='float64')
		#
		#
		#			# integrate each bin and interpolate to go from flux density to flux counts
		#			k_binspec.k_binspec(self.wave,self.atm,self.lambda_survey,Vatm)
		#			k_binspec.k_binspec(lam,		 sedphot, self.lambda_survey,Vsed)
		#			k_binspec.k_binspec(lam,		 sedflux, self.lambda_survey,Vflux)
		k_binspec.k_binspec(lam,		 flux, self.lambda_survey,flux_rebinned)
		#
		#
		#			# calculate mean
		#			median_Vsed  = numpy.median(Vsed)
		#			median_Vflux = numpy.median(Vflux)
		#
		#
		#			# find zeros
		#			check_Vsed = numpy.where(Vsed < 0.)[0]
		#			check_Vflux = numpy.where(Vflux < 0.)[0]
		#
		#
		#			# apply mean correction
		#			Vsed[check_Vsed] = median_Vsed
		#			Vflux[check_Vflux] = median_Vflux
		#				
		#
		#
		#			#=======================================================================================
		#			# Multiply fluxes by telescope transmission
		#			#   http://docs.scipy.org/doc/numpy/reference/generated/numpy.interp.html
		#			#=======================================================================================
		#			Vabsatm	= numpy.interp(self.lambda_survey,self.abswave,absatm)
		#			Vatm	= Vatm* Vtrans* (1.-Vabsatm)
		#			Vsed	= Vsed* Vtrans* (1.-Vabsatm)
		#
		#
		#
		#			#=======================================================================================
		#			# Calculate S/N
		#			#=======================================================================================
		#			Vtemp = Vatm+ Vsed+ self.num_readout* self.readnoise**2
		#			Vatmerr	 = numpy.sqrt(Vtemp)
		#			signoise = Vsed/Vatmerr
		#
		#			Vfluxerr = numpy.zeros(Vflux.shape)
		#			signoise_nonzero = numpy.where(signoise != 0)
		#			Vfluxerr[signoise_nonzero] = Vflux[signoise_nonzero]/signoise[signoise_nonzero]
		#				
		#			Vfluxsed = numpy.zeros(Vflux.shape)
		#			Vsed_nonzero = numpy.where(Vsed != 0)
		#			Vfluxsed[Vsed_nonzero] = Vflux[Vsed_nonzero]/Vsed[Vsed_nonzero]
		
		
		Vatmerr = numpy.sqrt((self.B + self.C*flux_rebinned)*(1 - self.D*air_mass_temp) + self.E)
		Vfluxerr = self.A*Vatmerr/(1-self.D*air_mass_temp)
		
		
		#=======================================================================================
		# Generate observed spectra
		#=======================================================================================
		random.seed(self.seed_constant)
		gasdev	  = numpy.random.exponential(scale=1,size=self.lambda_survey.size)
		Vfluxobs	= flux_rebinned + gasdev* Vfluxerr
		
		
		
		
		
		#			pylab.plot(self.lambda_survey,Vfluxerr,'r',self.lambda_survey,.999*Vfluxerr2,'g')
		#			pylab.show()
		
		return Vfluxobs, flux_rebinned	  #(Vflux is the the noise-less version)
		
		print Vfluxobs
	
	
	
	
	
	#=======================================================================================
	# Function to loop over galaxies
	#=======================================================================================
	def loop_over_gal(self):
		
		with archive.archive(data_bank,'r') as ar:
			# galaxy information
			self.flux_all			= ar['/gal/flux']
			self.lam_all			= ar['/gal/wavelength'] # one array only for every galaxy
			id_fiber			 = ar['/gal/fiber_id']
			fiber_selection_flag = ar['/gal/fiber_selection_flag']
		
		
		
		#=======================================================================================
		# cut data vectors
		#=======================================================================================
		where_id_fiber  = numpy.where(fiber_selection_flag == True)[0]
		id_fiber_temp   = id_fiber[where_id_fiber]
		self.air_mass	= self.air_mass[where_id_fiber]
		self.z			= self.z[where_id_fiber]
		nwhere_id_fiber = len(where_id_fiber)
		
		
		#=======================================================================================
		# LOOP over galaxies
		#=======================================================================================
		Ngal = len(id_fiber_temp)
		Vfluxobs_all = []
		Vfluxerr_all = []
		Vfluxobs_noisefree_all = []
		print "... Looping over", Ngal, "galaxies..."
		for i in range(Ngal):
			
			#print "... processing galaxy", i+1, "out of", Ngal
			flux 		  = numpy.array(self.flux_all[i], dtype='float32')
			lam			  = numpy.array(self.lam_all[i], dtype='float32')
			
			Vfluxobs, Vfluxobs_noise_free = self.observed_spectrum(lam, flux, i)
			
			# add to list
			Vfluxobs_all.append(Vfluxobs)
			Vfluxobs_noisefree_all.append(Vfluxobs_noise_free)
		
		print "... Looping over galaxies complete"
		
		
		
		#=======================================================================================
		# Output
		#=======================================================================================
		with archive.archive(data_bank,'a') as ar:
			ar['/gal/flux_spectrum_observed'] 			 = Vfluxobs_all
			ar['/gal/flux_spectrum_observed_noise_free'] = Vfluxobs_noisefree_all
			ar['/gal/wavelength_survey'] 				 = self.lambda_survey 		#!!! is this the right wavelength to save
		
		
		#=======================================================================================
		# DIAGNOSTIC
		# plot the last galaxy spectrum as an example
		#=======================================================================================
		
		print "... plotting diagnostics"
		
		try:
			Plotters.plot_spectrum(self.lambda_survey, Vfluxobs_noise_free,xmin=5000,xmax=10000,
								   fig_number=0, plot_number=0, title='Galaxy Spectrum Noise Free',
								   base_directory='.', filename_addendum="", show_plot=show_plot)
		except (RuntimeError, TypeError, NameError):
			print 'bad plot'
			pass
		
		try:
			Plotters.plot_spectrum(self.lambda_survey, Vfluxobs,xmin=5000,xmax=10000,
								   fig_number=1, plot_number=0, title='Galaxy Spectrum With Noise',
								   base_directory='.', filename_addendum="", show_plot=show_plot)
		except (RuntimeError, TypeError, NameError):
			print 'bad plot'
			pass




#=======================================================================================
# DIAGNOSTIC
#=======================================================================================
def main():
	
	s = SimulateObservedGalaxySpectra()
	s.loop_over_gal()
	print
	print

#-------------------------------
if __name__ == '__main__':
	main()
