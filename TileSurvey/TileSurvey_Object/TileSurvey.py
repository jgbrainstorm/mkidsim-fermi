import sys
sys.path.append('../../Wrapper/')
sys.path.append('../../Utilities/')
sys.path.append('../../Plotters/')
sys.path.append('../../Plotters_New/')

import archive
import Utilities
import Plotters
import Plotters_New

import random
import numpy



# =====================================================
# define throughput function
def survey_strategy_random():

	print
	print
	print "============================="
	print "		  Tiling Survey		 "
	print "============================="



	# =====================================================
	# Input
	# =====================================================
	databank_file = '../../../data/data_bank.h5'
	with archive.archive(databank_file,'r') as ar:
		# galaxies
		ra					  = ar['/gal/ra_true']
		dec					  = ar['/gal/dec_true']

		# tiles
		show_plot		 = ar['/Plotting/show_plot']

	seed = 9999

	print "... imported data"


	# count galaxies
	Ngal = len(target_selection_flag)


	# =====================================================
	# QA Test -- Incoming
	target_selection_flag_shape_expected = (Ngal,)
	check = (target_selection_flag.shape == target_selection_flag_shape_expected)
	if check :
		print "*** Incoming QA checksum passed for Survey Strategy***"
	else:
		print "XXX CHECKSUM FAILED XXX"
		print "target selection flag shape not compliant withe expectation"
		sys.exit(1)
	# =====================================================


	# set random seed
	random.seed(seed)


	# calculate number of galaxies to be in subset
	Nsub = int(Ngal*keep_fraction)


	# randomly sample a set of indices with the length of the subset
	sample = random.sample(numpy.arange(Ngal),Nsub)


	# create zeroed array
	survey_selection_flag = numpy.empty(Ngal,dtype=bool) 	#survey_selection_flag = numpy.zeros(Ngal)


	# set the sub sample to True
	survey_selection_flag[sample] = True


	# combine survey and target selection
	survey_selection_flag = survey_selection_flag* target_selection_flag
	print "... computed joint selection flag"

	# create list of tile_ID
	tile_id_tile = numpy.arange(Number_Tiles)

	# assign random values
	ramin = min(ra); ramax = max(ra)
	decmin = min(dec); decmax = max(dec)

	airmass = [random.uniform(1,1.01) for r in xrange(Number_Tiles)]
	tile_centers_ra = [random.uniform(ramin,ramax) for r in xrange(Number_Tiles)]
	tile_centers_dec = [random.uniform(decmin,decmax) for r in xrange(Number_Tiles)]
	#tile_centers_time = [random.uniform(timemin,timemax) for r in xrange(Number_Tiles)]


	# assign tile_ids
	check_survey_sel = numpy.where( survey_selection_flag  == True)[0]
	tile_id_galaxy = -numpy.ones(Ngal)
	tile_id_galaxy[check_survey_sel] = numpy.array([random.randint(0,Number_Tiles-1) for r in xrange(len(check_survey_sel))])
	print '... assigned tile ids'



	# =====================================================
	# QA Test --- Outgoing
	selected  = numpy.sum((survey_selection_flag == True))
	discarded = numpy.sum((survey_selection_flag == False))
	remainder = Ngal - selected - discarded
	if (remainder == 0):
		print "*** Outgoing QA checksum passed for Survey Strategy***"
	else:
		print "XXX CHECKSUM FAILED XXX"
		sys.exit(1)
	# =====================================================


	# count selected galaxies
	print '... gal counts'
	print '... ... n total',Ngal
	print '... ... n target selection',len( numpy.where( target_selection_flag == True)[0])
	print '... ... n both ', len( numpy.where( survey_selection_flag == True)[0])

	# compute min/max ra/dec
	ra_min_tot = numpy.min(ra[survey_selection_flag])
	ra_max_tot = numpy.max(ra[survey_selection_flag])
	dec_min_tot = numpy.min(dec[survey_selection_flag])
	dec_max_tot = numpy.max(dec[survey_selection_flag])



	# =====================================================
	# Output
	# =====================================================
	print '... output to databank'
	with archive.archive(databank_file,'a') as ar:
		# galaxies
		ar['/gal/survey_selection_flag'] = survey_selection_flag
		ar['/gal/tile_id']			   = tile_id_galaxy

		# tiles
		ar['/tiling/tile_id']		  = tile_id_tile
		ar['/tiling/number_of_tiles'] = Number_Tiles  # count tiles
		ar['/tiling/airmass']		  = airmass
		ar['/tiling/tile_centers_ra'] = tile_centers_ra
		ar['/tiling/tile_centers_dec'] = tile_centers_dec
		#ar['/tiling/tile_observation_time'] = tile_observation_time

		ar['/tiling/ra_min_tot'] = ra_min_tot
		ar['/tiling/ra_max_tot'] = ra_max_tot
		ar['/tiling/dec_min_tot'] = dec_min_tot
		ar['/tiling/dec_max_tot'] = dec_max_tot

	# =====================================================
	# Output
	# =====================================================
	print '... diagnostic plots'

	# ra
	try:
		width, center, bin_edges, hist = Utilities.make_histogram(tile_centers_ra, int(float(len(survey_selection_flag))/5.))
		xlabel = r"${\bf RA} {\rm [Deg]}$"
		ylabel = r"$N({\rm RA})$"
		Plotters_New.plot_line(center, hist, width=width,
					plot_type = 'both',
					#xmin = 0.,xmax=2.,
					xlabel=xlabel, ylabel=ylabel,
					color_bar = 'r',
					fig_name='ra_distribution',
					fig_number=0, plot_number=0, title='', base_directory='.', filename_addendum="", show_plot=show_plot
					)
	except (RuntimeError, TypeError, NameError):
		print 'bad plot; ra dec hist'
		pass

	# dec
	try:
		width, center, bin_edges, hist = Utilities.make_histogram(tile_centers_dec, int(float(len(survey_selection_flag))/5.))
		xlabel = r"${\bf DEC} {\rm [Deg]}$"
		ylabel = r"$N({\rm DEC})$"
		Plotters_New.plot_line(center, hist, width=width,
					plot_type = 'both',
					#xmin = 0.,xmax=2.,
					xlabel=xlabel, ylabel=ylabel,
					color_bar = 'r',
					fig_name='dec_distribution',
					fig_number=0, plot_number=0, title='', base_directory='.', filename_addendum="", show_plot=show_plot
					)
	except (RuntimeError, TypeError, NameError):
		print 'bad plot; ra dec hist'
		pass


	# airmass
	try:
		width, center, bin_edges, hist = Utilities.make_histogram(airmass, int(float(len(survey_selection_flag))/5.))
		xlabel = r"${\bf Airmass} $"
		ylabel = r"$N({\rm Airmass})$"
		Plotters_New.plot_line(center, hist, width=width,
					plot_type = 'both',
					xlabel=xlabel, ylabel=ylabel,
					color_bar = 'r',
					fig_name='airmass_distribution',
					fig_number=0, plot_number=0, title='', base_directory='.', filename_addendum="", show_plot=show_plot
					)
	except (RuntimeError, TypeError, NameError):
		print 'bad plot; ra dec hist'
		pass





#def importSurveyParameters(databank_file):

#	with archive.archive(databank_file,'r') as ar:
#		keep_fraction = ar['/Survey_Strategy/keep_fraction'] # fraction to keep in the dummy module
#		exposure_time = ar['/Survey_Strategy/exposure_time']
#		exposusre_delay = ar['/Survey_Strategy/exposure_delay']
#		time_after_dusk = ar['/Survey_Strategy/time_after_dusk']
#		num_bias_frame = ar['/Survey_Strategy/num_bias_frame']
#		num_flat_rate = ar['/Survey_Strategy/num_flat_rate']
#		location_alt = ar['/Survey_Strategy/location_alt']
#		location_lat = ar['/Survey_Strategy/location_lat']
#		location_lon = ar['/Survey_Strategy/location_lon']
#		airmass_max = ar['/Survey_Strategy/airmass_max']
#		mag_incr_back = ar['/Survey_Strategy/mag_incr_back']
#		seeing_mean = ar['/Survey_Strategy/seeing_mean']
#		seeing_sigma = ar['/Survey_Strategy/seeing_sigma']
#		ra_min = ar['/Survey_Strategy/ra_min']
#		ra_max = ar['/Survey_Strategy/ra_max']
#		dec_min = ar['/Survey_Strategy/dec_min']
#		dec_max = ar['/Survey_Strategy/dec_max']
#		date_start = ar['/Survey_Strategy/date_start']
#		date_end = ar['/Survey_Strategy/date_end']
#		include_moon_time = ar['/Survey_Strategy/include_moon_time']




# =====================================================
# Main
# =====================================================
def main():

	survey_strategy_random()
	print
	print


# =====================================================
if __name__ == '__main__':
	main()

