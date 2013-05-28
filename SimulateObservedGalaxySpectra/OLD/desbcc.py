import numpy as np
import scipy.interpolate as interpolate
import scipy.signal as signal
import scipy.mgrid as mgrid
import scipy.integrate as integrate
import archive
import random

# check the random to make sure it's gaussian random

def gauss_kern(sigma):
    """ Returns a normalized 1D gauss kernel array for convolutions """
    sigma = int(sigma)
    x = mgrid[-sigma:sigma+1]
    g = exp(-(x**2/float(sigma)))
    return g / g.sum()


def k_lambda_to_edges( lambda_centers):
	# N = num of pix centers
	# N+1 = num of pix edges

	nlambda				= len(lambda_centers)
	lambda_edges			= np.zeros(nlambda+1)
	lambda_edges[1:nlambda-1] 	= 0.5* (lambda_centers[0:nlambda-2]+ lambda_centers[1:nlambda-1L]) 
	lambda_edges[0]			= lambda_centers[0]- (lambda_edges[1]- lambda_centers[0])
	lambda_edges[nlambda]		= lambda_centers[nlambda- 1L]+(lambda_centers[nlambda- 1L]- lambda_edges[nlambda-1L])

	return lambda_edges



def k_smooth_py(loglam,flux,vdisp):

	nlambda		= len(loglam)
	pixsize		= np.abs(np.log(10.)*2.99792e5*(loglam[nlambda-1]-loglam[0])/float(nlambda))
	smoothing	= vdisp/pixsize    				# pixels
	npix		= long(4.0*np.ceil(smoothing))*long(2)+3
	klam		= findgen(npix)-float(npix-1.)/2.
	klam		= np.arange(float(npix))-float(npix-1.)/2.
	kernel		= exp(-0.5*(klam/smoothing)^2)/sqrt(2.*np.pi)/smoothing
	kernel		= kernel/sum(kernel)
	smoothed_spec	= signal.convol(flux,kernel)#,/edge_truncate)	 python code can't do this yet.


#def k_binspec_py(wave_in, spec_in, wave_survey):


def observe_galaxy():

	data_bank = 'data_bank.h5'
	with archive.archive(data_bank,'r') as ar:
		# galaxy information
		gal_air_mass    	= ar['/gal/airmass']

		# instrument information
		exposure 		= ar['/instrument/sensitivity/exposure_time']  # exposure = exposure time in seconds
		readnoise 		= ar['/instrument/sensitivity/read_noise']			# add to param files
		resnum 			= ar['/instrument/sensitivity/resolution']
		resolution_kms 		= ar['/instrument/sensitivity/resolution_kms']			# res_instrument = 300 # res of instrument  in km/s
		wavelength_min 		= ar['/instrument/wavelength_minimum']  	# Angstrom
		fiber_aperture_radius 	= ar['fiber_allocation/fiber_radius']
		worsening 		= ar['/instrument/sensitivity/worsening']			# add to param files
		altitude 		= ar['/instrument/sensitivity/altitude']			 #meters ??? add to param files

		# general constants
		plancks_constant 	= ar['/general/physical_constants/plancks_constant']	
		speed_of_light   	= ar['/general/physical_constants/speed_of_light']	

		# sky data
		wave   		= ar['sky/background/wave']
		temp    	= ar['sky/background/power']
		abswave 	= ar['sky/extinction/wave']
		abstemp 	= ar['sky/extinction/power']

		# Throughtput
		wavetrans 	= ar['/gal/throughput/wavelength']
		trans	  	= ar['/gal/throughput/throughput']  #!!! what should be the units here?


	
# !!! to modify to access tile information to link to where galaxy is observed
	#airmass=1.3

	Ngal = len(gal)

	Nwav=557*resnum		# 557 = number of pixels
	dispersion = 2.5		# def: dispersion = (pixsize or width ) / resolution; 2.5 is when most of the information comes 
					# in a limit
	dellam=7.14/resnum        #R is not being used. Only dellam matters now. = resolution/dispersion


	#=======================================================================================
	# Aperture
	#=======================================================================================
	aperture_area 	= np.pi*fiber_aperture_radius**2
	Area		= np.pi*4^2
	Areacm		= Area*10000.          #converting from m^2 to cm^2


	#=======================================================================================
	# SKY BACKGROUND
	# 	- Convert sky photons to photons/Ang.
	#	- in skybg_50_10.dat the units are photons/s/nm/arcsec^2/m^2
	#=======================================================================================
	atm	= worsening* exposure* temp* aperture_area* Area/10. #because 10ang=1 nm
	wave	= 10.* wave


	#=======================================================================================
	# TRANSMISSION information
	#=======================================================================================
	wavetrans = 10.* wavetrans		# WHAT are the translation units?
	trans	  = 0.01* trans

    
	#============================================
	# Modify ATMOSPHERIC ABSORPTION  
	#============================================
	absatm = abstemp*airmass*exp((1700.-altitude)/7000.)	
	absatm(np.where(absatm > 1)) = 1			 # file is bad ,...just make sure file is not bad
		

	#====================================
	# smooth spectra to SURVEY resolution
	# I think I shouldn't smooth over atmospheric absorption since 
	# smoothing only happens after the photons have been removed from atmosphere.
	#====================================
	logwave = np.log10(wave)
	atm 	= my_smooth(logwave, atm, resolution)  #atm 	= k_smooth(logwave,atm,300) #outatm=k_smooth(logwave,atm,1320)


	#=======================================================================================
	# Convert sed flux to photon counts.
	# 	- sed flux is in ergs/cm^2/s/Ang.
	# 	- photon counts is in photons/Ang.
	#=======================================================================================
	flux	= flux/ (1.+zspectemp) 			#add redshift dimming #flux[*,i]=flux[*,i]/(1.+zspectemp[i])^2 #add redshift dimming
	sedphot = exposure* flux* lam* Areacm/ (hp*c) 	#approximate - correct is to integrate 
	sedflux = exposure* flux* Areacm # 


	#=======================================================================================
	# Redshift galaxy sed's. 
	#	- rest-frame flux assumes galaxy is 10pc away.
	#	- doesn't matter when using coeffs calculated from
	#	- kcorrect, since the code already applies the 1/r^2 dimming.
	#=======================================================================================
	lamz = lam* (1. + zspectemp)

	#=======================================================================================
	# Interpolate wavelengths onto grid with resolution
	#	1) make grid 
	#	2) interpolate counts/fluxes onto grid
	#	3) Note that delta_lam will affect number of photons counted.
	#		***remember to multiply by delta_lam.***
	#=======================================================================================
	lambda_survey	 = np.zeros(Nwav)
	binwidth	 = np.zeros(Nwav)
	lambda_survey[0] = wavelength_minimum

	for i in np.arange(1,Nwave):
		lambda_survey[i] = lambda_survey[i-1]+dellam 

	for i in np.arange(1,Nwave-1):
		lambda_survey[i] = lambda_survey[i] + (lambda_survey[i+1]-lambda_survey[i])/2.
		binwidth[i] 	 = lambda_survey[i+1] - lambda_survey[i]

	#=======================================================================================
	#	Bin spectra to survey's resolutions by integrating over pixels
	#		to convert flux density to flux. 
	#=======================================================================================
	Vsed		= np.zeros(Nwav,Ngal)
	Vflux		= np.zeros(Nwav,Ngal)

	
	Vatm = k_binspec(wave,atm,lambda_survey)
	for i=0L,Nshort-1 do begin
		Vsed[*,i]=k_binspec(lamz[*,i],sedphot[*,i],lambda_survey)
		Vflux[*,i]=k_binspec(lamz[*,i],sedflux[*,i],lambda_survey)

	    
	#=======================================================================================
	##Multiply fluxes by telescope transmission
	#=======================================================================================
	Vtrans	= np.interp(lambda_survey,trans,wavetrans)
	Vabsatm	= np.interp(lambda_survey,absatm,abswave)
	Vatm	= Vatm* Vtrans* (1.-Vabsatm)
	Vsed	= Vsed* Vtrans* (1.-Vabsatm)


	#=======================================================================================
	##Calculate S/N
	#=======================================================================================
	Vatmerr	 = sqrt(Vatm+Vsed+readnoise**2)
	signoise = Vsed/Vatmerr
	Vfluxerr = Vflux/signoise


	#=======================================================================================
	##Generate observed spectra
	#=======================================================================================
	Vfluxobs=np.zeros(Nwav,Ngal)

	# loop over galaxies 
	for i in range(ngal):
		random.seed()
		gasdev 		= np.array([random.random() for x in range(Nwav)])
		Vfluxobs[*,i]	= Vflux[*,i]+ gasdev* Vfluxerr[*,i]

	# Make noise free
	Vfluxobs_noise_free=Vflux

	# Write out
	with archive.archive(data_bank,'a') as ar:
		ar['/gal/flux_spectrum_observed'] 	     = Vfluxobs
		ar['/gal/flux_spectrum_observed_noise_free'] = Vfluxobs_noise_free
		ar['/gal/wavelength_survey'] 		     = lambda_survey 		#!!! is this the right wavelength to save
		
