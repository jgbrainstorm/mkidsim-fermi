'''
(* applosses4.nb Will Saunders AAO Jan 2013. *)
(* 1) Sets up \
approximation for convolved profile for site seeing with Moffat \
\[Beta]=3.5 profile, exponential or deVaucouleurs R^1/4 image \
profiles, and Gaussian image degradation from spot sizes, guiding \
erros, dome/telescope seeing, deviations from ideal optics etc.
The approximation assumes a Moffat profile for the convolved profile, \
with \[Beta] and R chosen to give the correct central surface \
brightness and 2nd moment. *)
(* 2) Calculates aperture losses for \
this profile for circular apertures, including any offset from the \
image center. *)
(* 3) Estimates relative S/N at Zenith, assuming \
sky-limited observations, and including a factor for intrinsic \
spectral PSF degradation in the sperctrograph due to the finite fiber \
size. Plots and maximises S/N for various image types and parameters \
as a function of fiber size. *)
(* 4) calculates the differential \
refraction for arbitrary ZD and pair of wavelengths, according to ESO \
formula for LaSilla. Also the sky brightness and extinction as a \
fucnction of AM. *)
(* 5) Using SDSS AM distn, finds overall survey \
time to a fixed S/N for all target types,  tests these numbers for \
robsutness, and determines scalings between the input parameters and \
fiber size and survey time.*)
(* 6) Examines effect of ADC on survey \
time, including throughut derived from SK's ADC design.*)
(* To be read in conjunction with technical note 'DESpec_ADC \
_aperture _survey _speed', WS 1/2013 *)

'''

import sys
sys.path.append('../../../Wrapper/')
import archive
import pylab
import numpy
from scipy.optimize import root, fsolve
from scipy.integrate import quad, tplquad





# ============================================
# Utilities
# ============================================
def check_match(check_answer, solution, tolerance=1e-4, check_type="", func_name=""):
	diff = abs(solution - check_answer)
	diff_frac = diff/solution
	if diff_frac >= tolerance : print "uh oh, no match ("+func_name+", "+check_type+")", solution, check_answer

	return

def check_profiles():
	print '... checking profiles'
	Moff35h_height_check(func_name='moff35h height', check_type='height')
	integrate_Moff(func_name='moff', check_type='integrate')
	integrate_MoffM2(func_name='moffm2', check_type='integrate')
	integrate_check_gauss(func_name='gauss', check_type='integrate')
	integrate_check_ae(func_name='ae', check_type='integrate')
	integrate_check_eprof(func_name='eprof', check_type='integrate')
	integrate_check_vproft(func_name='vproft', check_type='integrate')
	integrate_check_vprof(func_name='vprof', check_type='integrate')


# ============================================
# Moffat Profile
# ============================================
# --- Basic Moffat Profile
def Moff(r, R, beta):
	beta_a = beta - 1.
	rr	   = 1. + (r/R)**2
	bottom = numpy.pi* R**2* rr**beta
	result = beta_a / bottom

	return result


# --- Integrate Moffat to check normalization
def integrate_Moff(R=1.2345, beta=3.5567, check_type=None, func_name=None ):
	 moffat_psf_integrated = quad(lambda r: 2.* numpy.pi* r* Moff(r, 1.2345, 3.5567), 0, numpy.inf)[0]
	 check_match(1., moffat_psf_integrated,check_type=check_type,func_name=func_name)

	 return


# --- Second moment of Moffat profile
def MoffM2(r, RM, beta):
	beta_a = beta- 1.
	beta_b = beta- 2.
	rrM	   = 1.+ ((r/ RM)**2/ beta_b)
	bottom = numpy.pi* beta_b* RM**2* rrM**beta
	result = beta_a / bottom

	return result


# --- Integrate Moffat second moment to check normalization
# --- ---  so for MoffM2, M2 = RM^2 always
def integrate_MoffM2(R=1.2345, beta=3.5678, check_type=None, func_name=None ):
	moffat_psf_2mom_integrated = quad(lambda r: 2.* numpy.pi* r* r**2* MoffM2(r, 1.2345, 3.5567)/ R**2, 0, numpy.inf)[0]
	check_match(1., moffat_psf_2mom_integrated,check_type=check_type,func_name=func_name)
	return

# --- find HWHM
# --- --- r /. FindRoot[Moff[r, 1., 3.5] == Moff[0., 1., 3.5]/2, {r, .5}]
def Moffroot(r):
	result = Moff(r,R=1.,beta=3.5) - Moff(0., R=1., beta=3.5)/2.

	return result


def Rhwhm35(rguess = 0.5):
	check_answer = 0.467989
	solution = fsolve(Moffroot,rguess, xtol=1e-06, maxfev=5000)[0]
	check_match(check_answer, solution)
	#print 'check rhwhm35', solution

	return solution


# --- Moffat with be = 3.5
# --- ---  Moff35h is a Moffat profile with beta=3.5, defined by its half-width at half height
def Moff35h(r, Rhwhm, Rhwhm35=0.467989) :
	fac_a = 1./(2.* numpy.pi) * (Rhwhm**2 / Rhwhm35**2)
	fac_b = ( 1. + (Rhwhm35* r/ Rhwhm)**2 )**(-3.5)
	result  = 5.* fac_a* fac_b

	return result


# --- check half-height
def Moff35h_height_check(check_type=None, func_name=None):
	result = Moff35h(1, 1)/ Moff35h(0, 1)
	check_match(0.5, result,check_type=check_type,func_name=func_name)

	return


# --- --- FindRoot[Integrate[2 Pi r r^2 Moff35h[r, 1], {r, 0, Ref}] == 0.5, {Ref, 1}]
def integrand_moff(Ref, Compare = 0.5):
	result = 2.* numpy.pi* quad(lambda r : r* r**2* Moff35h(r, 1.), 0., Ref, epsabs=1.e-9,epsrel=1.e-9)[0] - Compare
	return result

def root_check_moff(Ref_guess=1.):
	solution = fsolve(integrand_moff,Ref_guess, xtol=1e-09, maxfev=9000)[0]
	#print '***root moff***',solution, "(should be 1.44585)" # !!!
	return solution

# ============================================
# Gaussian Profile
# ============================================
def GProf(r, Rrms):
	bottom = numpy.pi* Rrms**2
	top	   = numpy.exp(-(r/Rrms)**2)
	result = top / bottom

	return result


def integrate_check_gauss(check_type=None, func_name=None):
	gauss_integrated = quad(lambda r: 2.* numpy.pi* r* GProf(r,1.2345), 0, numpy.inf)[0]
	check_match(1., gauss_integrated,check_type=check_type,func_name=func_name)

	return


# ============================================
# Exponential Profile
#   ... where to get ae
# ============================================
def EProf(r, Re50, ae=1.67835):
	fac_a = 1./(numpy.pi*2.)
	fac_b = (ae/Re50)**2
	fac_c = numpy.exp(-ae*r/Re50)
	result = fac_a* fac_b* fac_c

	return result


# --- --- ae = r /. FindRoot[(1 + r) Exp[-r] == 0.5, {r, 1}]
def ae(rguess = 1., Compare = 0.5):
	ae = fsolve(lambda r: (1.+ r)* numpy.exp(-r) - Compare,rguess, xtol=1e-06, maxfev=5000)[0]
	check_match(1.67835, ae)

	return ae


def integrate_check_ae(check_type=None, func_name=None):
	ae_integrated = quad(lambda r: ae()**2* r* numpy.exp(-ae()* r), 0., 1.)[0]
	check_match(0.5, ae_integrated,check_type=check_type,func_name=func_name)

	return ae_integrated


# --- --- Integrate[2 Pi r EProf[r, 1], {r, 0, numpy.infinity}]
def integrate_check_eprof(check_type=None, func_name=None):
	eprof_integrated = quad(lambda r: 2.* numpy.pi* r* EProf(r,1), 0, numpy.inf)[0]
	check_match(1., eprof_integrated,check_type=check_type,func_name=func_name)

	return eprof_integrated




# ============================================
# deVaucoleurs Profile
# ============================================
def VProft_original(r, Re50):
	fac_a = (1./ Re50)**2
	fac_b = numpy.exp(-(r/Re50)**0.25)
	result = fac_a* fac_b

	return result

def VProft(r, Re50, KV=126669.):
	fac_a = (1./ Re50)**2/ KV
	fac_b = numpy.exp(-(r/Re50)**0.25)
	result = fac_a* fac_b

	return result

def VProf(r, Re50, Reft=3459.4850148369196, KV=126669.):
	fac_a = (Reft/ Re50)**2/ KV
	fac_b = numpy.exp(-(Reft* r/Re50)**0.25)
	result = fac_a* fac_b

	return result


# --- --- KV = Integrate[2 Pi r VProft[r, Ref], {r, 0, numpy.infinity}, Assumptions -> Re[Ref^0.25] > 0]
def KV():
	Ref = 1.
	result = quad(lambda r : 2.* numpy.pi* r* VProft_original(r, Ref), 0, numpy.inf)[0]
	check_match(126669., result)
	return result

def integrate_check_vproft(check_type=None, func_name=None):
	Ref = 1.
	result = quad(lambda r : 2.* numpy.pi* r* VProft(r, Ref), 0, numpy.inf)[0]
	check_match(1., result,check_type=check_type,func_name=func_name)
	return


# --- --- Reft = Ref /. FindRoot[Integrate[2 Pi r  VProft[r, 1], {r, 0, Ref}] == 0.5, {Ref,1}]
def Integrand_Reft(Ref, Compare=0.5):
	Ref_constant = 1.
	result = quad(lambda r : 2.* numpy.pi* r* VProft(r, Ref_constant), 0, Ref)[0] - Compare
	return result

def Reft(rguess=1):
	solution = fsolve(Integrand_Reft,rguess, xtol=1e-06, maxfev=5000)[0]
	return solution

def integrate_check_vprof(check_type=None, func_name=None):
	Ref = 1.
	result = quad(lambda r : 2.* numpy.pi* r* VProf(r, Ref), 0, numpy.inf)[0]
	check_match(1., result,check_type=check_type,func_name=func_name)
	return



# ============================================
# Second order moments
# ============================================

# --- --- (* So the 2nd moment of Moff35[r,Rhwhm] is 3.044 Rhwhm^2 *)
def MM2():  # Moffat
	print '... ... ... mm2'
	upper_limit = numpy.inf
	#Rnew = 1.  # !!!
	Rnew = 0.573165 # !!!
	result = quad(lambda r: 2.* numpy.pi* r* r**2* Moff35h(r, Rnew), 0., upper_limit, epsrel=1e-9, epsabs=1e-9,limit=1000, limlst=1000)[0]
	#check_match(3.04395, result, func_name='mm2') !!!

	return result


def MG2():  # Gauss
	print '... ... ... mg2'
	result = quad(lambda r: 2.* numpy.pi* r* r**2* GProf(r, 1), 0., numpy.inf)[0]
	check_match(1., result, func_name='mg2')

	return result


def MV2():  # deVaucoleurs
	print '... ... ... mv2'
	result = quad(lambda r: 2.* numpy.pi* r* r**2* VProf(r, 1), 0., numpy.inf)[0]
	check_match(21.6793, result, func_name='mv2')

	return result


def ME2():  # Exponential
	print '... ... ... me2'
	result = quad(lambda r: 2.* numpy.pi* r* r**2* EProf(r, 1), 0., numpy.inf)[0]
	check_match(2.13004, result, func_name='me2')

	return result

#!!!
def MMEVG2(Rhwhm, RefE, RefV, Rrms) :
	print '... ... mmevg2'
	result = MM2()* Rhwhm**2 + ME2()* RefE**2 + MV2()* RefV**2 + MG2()* Rrms**2
	#result = MM2()* Rhwhm**2
	#result = ME2()* RefE**2
	#result = MV2()* RefV**2
	#result = MG2()* Rrms**2
	return result


# ============================================
# Plot Profiles
# ============================================
def plot_profiles(fig_number=0, plot_number=0,  title='apertures'):
	print "... plot profiles"

	fig = pylab.figure(fig_number)
	file_figure_out = "./Figures/aperture_profiles_"+str(plot_number)+".png"

	pylab.xlabel('radius ')
	pylab.ylabel('height')
	pylab.title(title)
	pylab.xticks(fontsize=10)
	pylab.yticks(fontsize=10)
	pylab.grid(True)

	radius = numpy.arange(0.,4,0.01)
	moffat = Moff(radius, 1., 3.5)
	moff35h = Moff35h(radius, 1.)
	gauss  = GProf(radius, 1.)
	expon  = EProf(radius, 1., ae=1.)
	devauc = VProf(radius, 1.)

	pylab.plot(radius, moffat, 'r-', label='moffat')
	pylab.plot(radius, moff35h, 'c-', label='moff35h')
	pylab.plot(radius, gauss , 'b-', label='gauss')
	pylab.plot(radius, expon , 'k-', label='expon')
	pylab.plot(radius, devauc, 'g-', label='devauc')
	#pylab.loglog()

	pylab.legend(loc=0,prop={'size':10})
	pylab.savefig(file_figure_out)
	print "... ... plot saved to ", file_figure_out



# ============================================
# MG0
# ============================================
def MG0(Rhwhm, Rrms):
	result = quad(lambda r: 2.* numpy.pi* r* MoffM2(r, 1.7446920204827057* Rhwhm , 3.5)* GProf(r,Rrms),
					0, 3*Rrms,
					epsrel=3. )[0]

	return result


# ============================================
# MEVG0 is the value of the MEG or MVG convolution at the origin
# ============================================
def MEVG0_integrand_exp(r, rt, th, Rhwhm, RefE, RefV, Rrms):
	print '... ... ... mevg0 integrand'
	result = \
		2.* numpy.pi* r*\
		MoffM2(numpy.sqrt(r**2+ rt**2- 2* r* rt* numpy.cos(th)), 1.7446920204827057* Rhwhm, 3.5)*\
		2* rt*\
		GProf(rt, Rrms)*\
		EProf(r,RefE)
	return result

def MEVG0_integrand_vau(r, rt, th, Rhwhm, RefE, RefV, Rrms):
	result = \
		2.* numpy.pi* r*\
		MoffM2(numpy.sqrt(r**2+ rt**2- 2* r* rt* numpy.cos(th)), 1.7446920204827057* Rhwhm, 3.5)*\
		2* rt*\
		GProf(rt, Rrms)*\
		VProf(r, RefV),
	return result



def MEVG0(Rhwhm, RefE, RefV, Rrms) :
	print '... ... ... MEVG0: Rhwhm, RefE, RefV, Rrms', Rhwhm, RefE, RefV, Rrms

	if RefE !=0 and RefV !=0 :
		print "Ooops RE OR RV please"
		sys.exit(0)

	if RefE == 0 and RefV == 0:
		result = MG0(Rhwhm, Rrms)

	if RefV == 0 and RefE != 0 :
		result = tplquad(MEVG0_integrand_exp,\
						0.0, numpy.pi, \
						lambda th: 0.0,\
						lambda th: 3.0*Rrms,\
						lambda th,rt: 0.0,\
						lambda th,rt: 5.0*RefE,\
						args = (Rhwhm, RefE, RefV, Rrms),\
						)[0]

	if RefE == 0 and RefV != 0:
		result = tplquad(MEVG0_integrand_vau,\
						0, numpy.pi,
						lambda th: 0.0,\
						lambda th: 3.0*Rrms,\
						lambda th,rt: 0.0,\
						lambda th,rt: 5.0*RefV,\
						args = (Rhwhm, RefE, RefV, Rrms),\
						)[0]


	return result




# ============================================
# ============================================
def BetaMEVG0(Rhwhm, RefE, RefV,Rrms) :
		print "... ... betamevg0"
		fac_a = (2*numpy.pi* MEVG0(Rhwhm, RefE, RefV, Rrms)* MMEVG2(Rhwhm, RefE, RefV, Rrms) - 1)
		fac_b = (  numpy.pi* MEVG0(Rhwhm, RefE, RefV, Rrms)* MMEVG2(Rhwhm, RefE, RefV, Rrms) - 1)
		result = fac_a / fac_b

		return result


# ============================================
# constants that go into ApEtaAll
# ============================================
# Environmental Constants for La Silla
TC = 11.5		  # temperature assumptino [celsius]
T  = 273.16 + TC   # temperature Kelvin
PS = -10474.0 + 116.43* T - 0.43284* T**2 + 0.00053840* T**3
RH = 0.442  # relative humidity
P2 = RH* PS # pressure
P1 = 772.2  # pressure

# factors for galaxies
Relg	= 0.4
Relgmax = 0.5
Rlrg	= 0.35
Rlrgmax = 0.5

# Plate scale (despec)
PlateScDESP = 56.73			 # micron/arcsec
FWHMCT		= 0.75				# seeing at cTIO
FWHMbyrms   = 2.* numpy.sqrt(numpy.log(4.))/numpy.sqrt(2.)
FWHMfixDS   = 0.5



# ============================================
# functions #???
# ============================================
def D1(P1, T):
	result = P1/T*\
			(1.0 +\
				P1*\
				(57.90*  1e-88 - (9.3250* 1e-4/ T) + (0.25844/T**2))\
			)
	return result

def D2(P2, T):
	result = P2/T*\
			(1.0 +\
				P2*\
				(1.0+ 3.7*1e-4* P2)*\
				(-2.37321* 1e-3 + (2.23366/ T) - (710.792/T**2) + (7.75141* 1e4/T**3))\
			)
	return result

def N1(wave):
	fac_a = (2371.34 + 683939.7/ (130.- 1./ wave**2) + 4547.3/(38.9 - 1./wave**2))
	fac_b = (6487.31 + 58.058/wave**2 - 0.71150/wave**4 +   0.08851/wave**6)
	result = 1e-8 * fac_a* D1(P1,T) + fac_b* D2(P2,T)

	return result


def ZD(AM): # zenith distance
	return numpy.arccos(1/AM)


def DR(wave, wave_0, AM):
	fac_a = numpy.tan(ZD(AM))
	fac_b = (N1(wave_0/1000.) - N1(wave/1000.))* 206264.8
	result = fac_a * fac_b

	return result


def Rrmsfix(FWHMfixDS, PlateScDESP, FWHMbyrms):
	result = FWHMfixDS* PlateScDESP/ FWHMbyrms

	return result


def Rsee(FWHMsite, wave, AM):
	fac_a = FWHMsite/2.
	fac_b = (wave/500.)**(-0.2)
	fac_c = AM**(3./5.)
	result = fac_a* fac_b* fac_c

	return result

# ============================================
# Apetad (defocused)
# ============================================
def Frac(r1, r2, delta=0) :
	if delta == 0 :
		result = 0
	else :
		dist0 = (r1**2 + delta**2 - r2**2 )/ (2*r1*abs(delta))
		min0  = min(1,dist0)
		max0  = max(-1,min0)
		result = (numpy.arccos(max0)/ numpy.pi)
		result = result.real

	return result


def ApEtad(R, Beta, RAp, delta=0) :
	print '... ... apetad'
	if delta == 0:
		result = quad(lambda r: 2* numpy.pi* r* MoffM2(r, R, Beta),  0, RAp)[0]
	else :
		result = quad(lambda r: 2* numpy.pi* r* Frac(r, RAp, delta)* MoffM2(r, R, Beta),  0, numpy.inf)[0]

	return result


# ============================================
# ApEtaAll
# ============================================
'''
 ApEtaAll is general purpose aperture loss function
		 \[Lambda],\[Lambda]0 are target wavelength and fiber-centering \
				 wavelength
				 AM is airmass
				 Decent is decentering in um (added in quadrature to differential \
						 refraction)
				 RmsSpot is spot size from ZEMAX, assumed Gaussian
				  FWHMFix is addition for dome seeing, guiding, alignment etc
				  DDefoc is defocus diameter = axial dist/f-number
				  PlateSc is plate scale (um/arcsec)
				  DAp is fiber diamater
				  RE is effective radius for exponential profile
				  RV is effective radius for deVaucouleurs profile
'''
def ApEtaAll(wave, wave0, FWHMsite, AM, Decent, RmsSpot, FWHMFix, DDefoc, PlateSc, DAp, RE, RV) :
	print 'starting apetaall'
	print '... fwhm'
	FWHM_temp  = numpy.sqrt((FWHMFix/FWHMbyrms)**2 + (RmsSpot/PlateSc)**2 + (DDefoc/PlateSc)**2/8)
	print '... ...', FWHM_temp

	print '... Rsee'
	Rsee_temp  = Rsee(FWHMsite, wave, AM)
	print '... ...', Rsee_temp

	print '... d-temp'
	d_temp	 = numpy.sqrt((Decent/PlateSc)**2 + DR(wave, wave0, AM)**2)
	print '... ...', d_temp

	print '... MMEVG'
	MMEVG_temp = numpy.sqrt(MMEVG2(Rsee_temp, RE, RV, FWHM_temp))
	print '... ...', MMEVG_temp

	print '... betamevg0'
	BetaMEVG0_temp =  BetaMEVG0(Rsee_temp, RE, RV, FWHM_temp)
	print '... ...', BetaMEVG0_temp

	print '... apetad'
	result = ApEtad( MMEVG_temp, BetaMEVG0_temp, DAp/2, d_temp )
	print '... ...', result
	print
	sys.exit(1)

	return result


# ============================================
# Main
# ============================================
def main():
    # wave      :   target wavelength [nm]
    # wave 0    :   fiber-centering wavelength [nm]
    # AM        :   airmass
    # Decent    :   decentering [um]
    # RmsSpot   :   spot size (from ZEMAX, assumed Gaussian)
    # FWHMFix   :   addition to dome seeing, etc.
    # DDefoc    :   defocus diameter = axial distance/f-number
    # PlateSc   : plate scale [um/arcsec]
    # DAp       :   fiber diameter [mm] # ???
    # RE        :   effective radius for exponential profile
    # RV        :   effective radius for deVaucouleurs profile

	plot_profiles()
	check_profiles()

	wave, wave0, FWHMsite, AM, Decent, RmsSpot, FWHMFix, DDefoc, PlateSc,	DAp, RE, RV = \
	350, 400,	FWHMCT,  1.0, 0,	  25,	  0.5,		20, PlateScDESP, 1.76, 0, 0
		#350, 400, FWHMCT, 1.0, 10, 25, 0.5, 20, PlateScDESP, 1.76, 0, 0

	result = ApEtaAll(wave, wave0, FWHMsite, AM, Decent, RmsSpot, FWHMFix, DDefoc, PlateSc, DAp, RE, RV)


#-------------------------------
if __name__ == '__main__':
		main()
