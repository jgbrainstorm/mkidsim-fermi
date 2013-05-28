# Throughput for DESpec
# Created by Brian Nord @ FNAL, nord@fnal.gov, July, 2012
# Derived from work by Will Saunders

# import
import sys
sys.path.append('../../Wrapper/')
sys.path.append('../../Plotters/')

import archive
import Plotters

import numpy


# =====================================================
# define throughput function
def throughput():

	# =====================================================
	# INPUT
	# =====================================================
	databank_file = '../../../data/data_bank.h5'

	print
	print
	print "========================="
	print "entering Throughput"
	print "========================="


	with archive.archive(databank_file, 'r') as ar:
		# --------------------------------------------------------
		# read hardcoded through puts for large elements
		wavelength_temp			= ar['/Throughput/parameters/wavelength']
		aperture_losses_gal 	= ar['/Throughput/parameters/aperture_losses_gal']
		aperture_losses_star 	= ar['/Throughput/parameters/aperture_losses_star']
		mohawk_frd  			= ar['/Throughput/parameters/mohawk_frd']
		collimator 				= ar['/Throughput/parameters/collimator']
		vph_gsolver 			= ar['/Throughput/parameters/vph_gsolver']
		camera 					= ar['/Throughput/parameters/camera']
		ccd						= ar['/Throughput/parameters/ccd']

		# --------------------------------------------------------
		# read hardcoded throuhg puts from glass/coatings
		Al 			= ar['/Throughput/parameters/Al']
		B270_air_glass 	= ar['/Throughput/parameters/B270_air_glass']
		air_Bk7 	= ar['/Throughput/parameters/air_Bk7']
		B270_6mm 	= ar['/Throughput/parameters/B270_6mm']
		ctio 		= ar['/Throughput/parameters/ctio']
		air_silica 	= ar['/Throughput/parameters/air_silica']
		blue_MgF 	= ar['/Throughput/parameters/blue_MgF']
		fiber_material 	= ar['/Throughput/parameters/fiber_material']
		red_MgF 	= ar['/Throughput/parameters/red_MgF']
		LLF1 		= ar['/Throughput/parameters/LLF1']
		seso_a 		= ar['/Throughput/parameters/seso_a']
		seso_b 		= ar['/Throughput/parameters/seso_b']
		sf5 		= ar['/Throughput/parameters/sf5']
		bk7_25mm 	= ar['/Throughput/parameters/bk7_25mm']
		lf5 		= ar['/Throughput/parameters/lf5']
		fk5 		= ar['/Throughput/parameters/fk5']
		lak33 		= ar['/Throughput/parameters/lak33']
		laf21 		= ar['/Throughput/parameters/laf21']
		prot_ag 	= ar['/Throughput/parameters/prot_ag']
		ag_a		= ar['/Throughput/parameters/ag_a']
		ag_b 		= ar['/Throughput/parameters/ag_b']
		blue_broad 	= ar['/Throughput/parameters/blue_broad']
		red_broad 	= ar['/Throughput/parameters/red_broad']
		edmond 		= ar['/Throughput/parameters/edmond']
		solgel_b 	= ar['/Throughput/parameters/solgel_b']
		solgel_r 	= ar['/Throughput/parameters/solgel_r']
		solgel_plus	= ar['/Throughput/parameters/solgel_plus']

	print "... imported data"

	# --------------------------------------------------------
	# defining each optical element in terms of material
	# throughput (label with column in original worksheet) ... can replace the some of the above
	#atmos 		= 10**(-0.4*ctio*1.15)
	primary 	= 0.98 * Al
	top_end 	= 1.-(2.91/(numpy.pi/4.* 4.**2))   	#
	wfc			= air_silica**2* seso_b**6* red_broad**4* 10**(-fiber_material* 0.000376 /10.)
	#adc 		= red_broad**4 *LLF1**(28/25)*bk7_25mm**(40.5/25)	#
	fiber 		= 10**(-fiber_material*0.05/10)				#
	print "... Combined minor optical elements"


	# --------------------------------------------------------
	# combining optical elements
	telescope 	 = primary * top_end * wfc
	mohawk_full	 = fiber * mohawk_frd			#
	vph_mounted	 = vph_gsolver * red_broad**2 * bk7_25mm**(20./25.)
	spectrograph = collimator * vph_mounted * camera * ccd
	print "... Combined major optical elements with galaxy, including aperture losses"

	# --------------------------------------------------------
	# combine elements with aperture losses for gal and star
	# this will have to be in a loop when applying the aperture losses for each gal
	#		  for a star, total_star 	= telescope * aperture_losses_star * mohawk_full * spectrograph
	#	will have to apply galaxy aperture in SIMGALSPECOBSERVED
	total_gal 	= telescope * aperture_losses_gal  * mohawk_full * spectrograph
	throughput_final = total_gal
	print "... computed final throughput"


	# =====================================================
	# Output
	# =====================================================
	with archive.archive(databank_file, 'a') as ar:
		ar['/Throughput/throughput'] = throughput_final

	print "... exported data to databank"

	# =====================================================
	# Diagnostic
	# =====================================================
	print "... basic diagnostic plots"
	Plotters.plot_throughput_spectrum(wavelength_temp, throughput_final, power1=telescope, power2=spectrograph, xmin=500, xmax=1000,
							fig_number=0, plot_number=0,
							label0 = '', label1 = '', label2='',
							title='ThroughputSpectrum',
							base_directory='.', filename_addendum="", show_plot=False)

	return

# =====================================================
# Main
# =====================================================
def main():

	throughput()
	print
	print


#-------------------------------
if __name__ == '__main__':
	main()


