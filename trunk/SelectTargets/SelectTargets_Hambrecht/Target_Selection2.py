import sys
sys.path.append('../../Wrapper/')
import archive
sys.path.append('../../Plotters/')
import Plotters
import random
import numpy
import scipy.interpolate
import csvimport
from matplotlib import pyplot
import os

print
print
print "============================="
print "entered Target Selection"
print "============================="


inputfile = '../../../data/data_bank.h5'
outputfile = inputfile

# import parameters
with archive.archive(inputfile,'r') as ar:
	training_catalog_location = ar['/Target_Selection/training_catalog_location']
	mag_g_min = ar['/Target_Selection/mag_g_min']
	mag_g_max = ar['/Target_Selection/mag_g_max']
	mag_r_min = ar['/Target_Selection/mag_r_min']
	mag_r_max = ar['/Target_Selection/mag_r_max']
	mag_i_min = ar['/Target_Selection/mag_i_min']
	mag_i_max = ar['/Target_Selection/mag_i_max']
	mag_z_min = ar['/Target_Selection/mag_z_min']
	mag_z_max = ar['/Target_Selection/mag_z_max']
	mag_y_min = ar['/Target_Selection/mag_y_min']
	mag_y_max = ar['/Target_Selection/mag_y_max']
	photo_z_min = ar['/Target_Selection/photo_z_min']
	photo_z_max = ar['/Target_Selection/photo_z_max']
	type_cutoff = ar['/Target_Selection/type_cutoff']
	max_training_set_size = ar['/Target_Selection/max_training_set_size']
	recalculate_weights = ar['/Target_Selection/recalculate_weights']
	weightsfile = ar['/Target_Selection/weights_file_location']

if (not recalculate_weights):
	with archive.archive(weightsfile,'r') as ar:
		weights = ar['/weights']

# import mags and photo-zs from target cat
with archive.archive(inputfile,'r') as ar:
	magnitude_g = ar['/gal/magnitude_g']
	magnitude_r = ar['/gal/magnitude_r']
	magnitude_i = ar['/gal/magnitude_i']
	magnitude_z = ar['/gal/magnitude_z']
	magnitude_y = ar['/gal/magnitude_y']
	photo_z = ar['/gal/redshift_photo_gaussian']

print "... finished data param import"


numberOfGalaxies = magnitude_y.shape[0]




if (training_catalog_location == ''):
	# we'll use our VISTA cat and our selection parameters

	# cut out the 99s and -99s
	magcut99_g_flag = (magnitude_g < 99)*(magnitude_g > -99)
	magcut99_r_flag = (magnitude_r < 99)*(magnitude_r > -99)
	magcut99_i_flag = (magnitude_i < 99)*(magnitude_i > -99)
	magcut99_z_flag = (magnitude_z < 99)*(magnitude_z > -99)
	magcut99_y_flag = (magnitude_y < 99)*(magnitude_y > -99)
	magcut99_flag = magcut99_g_flag*magcut99_r_flag*magcut99_i_flag*magcut99_z_flag*magcut99_y_flag

	magnitude_g_99cut = magnitude_g[magcut99_flag]
	magnitude_r_99cut = magnitude_r[magcut99_flag]
	magnitude_i_99cut = magnitude_i[magcut99_flag]
	magnitude_z_99cut = magnitude_z[magcut99_flag]
	magnitude_y_99cut = magnitude_y[magcut99_flag]
	photo_z_99cut = photo_z[magcut99_flag]


	# perform mag and photo-z cuts
	magcut_g_flag = (magnitude_g > mag_g_min)*(magnitude_g < mag_g_max)
	magcut_r_flag = (magnitude_r > mag_r_min)*(magnitude_r < mag_r_max)
	magcut_i_flag = (magnitude_i > mag_i_min)*(magnitude_i < mag_i_max)
	magcut_z_flag = (magnitude_z > mag_z_min)*(magnitude_z < mag_z_max)
	magcut_y_flag = (magnitude_y > mag_y_min)*(magnitude_y < mag_y_max)

	magcut_flag = magcut_g_flag*magcut_r_flag*magcut_i_flag*magcut_z_flag*magcut_y_flag

	photozcut_flag = (photo_z > photo_z_min)*(photo_z < photo_z_max)

	simplecuts_flag = magcut_flag*photozcut_flag

	magnitude_g_cut = magnitude_g[simplecuts_flag]
	magnitude_r_cut = magnitude_r[simplecuts_flag]
	magnitude_i_cut = magnitude_i[simplecuts_flag]
	magnitude_z_cut = magnitude_z[simplecuts_flag]
	magnitude_y_cut = magnitude_y[simplecuts_flag]
	photo_z_cut = photo_z[simplecuts_flag]


	print "... performed simple cuts"

	if (recalculate_weights):
		print "... Training..."

		# import from training catalog
		with archive.archive(inputfile,'r') as ar:
			magnitude_g_training = ar['/target_selection_training/gal/magnitude_g']
			magnitude_r_training = ar['/target_selection_training/gal/magnitude_r']
			magnitude_i_training = ar['/target_selection_training/gal/magnitude_i']
			magnitude_z_training = ar['/target_selection_training/gal/magnitude_z']
			magnitude_y_training = ar['/target_selection_training/gal/magnitude_y']
			photo_z_training	 = ar['/target_selection_training/gal/redshift_photometric']
			galaxy_type_training = ar['/target_selection_training/gal/galaxy_type']

		print "... ... finished import of training cat"


		# perform mag and photo-z cuts
		magcut_g_training_flag = (magnitude_g_training > mag_g_min)*(magnitude_g_training < mag_g_max)
		magcut_r_training_flag = (magnitude_r_training > mag_r_min)*(magnitude_r_training < mag_r_max)
		magcut_i_training_flag = (magnitude_i_training > mag_i_min)*(magnitude_i_training < mag_i_max)
		magcut_z_training_flag = (magnitude_z_training > mag_z_min)*(magnitude_z_training < mag_z_max)
		magcut_y_training_flag = (magnitude_y_training > mag_y_min)*(magnitude_y_training < mag_y_max)

		magcut_training_flag = magcut_g_training_flag*magcut_r_training_flag*magcut_i_training_flag*magcut_z_training_flag*magcut_y_training_flag

		photozcut_training_flag = (photo_z_training > photo_z_min)*(photo_z_training < photo_z_max)

		simplecuts_training_flag = magcut_training_flag*photozcut_training_flag

		magnitude_g_training = magnitude_g_training[simplecuts_training_flag]
		magnitude_r_training = magnitude_r_training[simplecuts_training_flag]
		magnitude_i_training = magnitude_i_training[simplecuts_training_flag]
		magnitude_z_training = magnitude_z_training[simplecuts_training_flag]
		magnitude_y_training = magnitude_y_training[simplecuts_training_flag]
		photo_z_training = photo_z_training[simplecuts_training_flag]
		galaxy_type_training = galaxy_type_training[simplecuts_training_flag]


		print "... ...  performed simple cuts on training cat"

		#   cap to max_training_set_size
		if (numpy.sum(simplecuts_training_flag) > max_training_set_size):
			magnitude_g_training = magnitude_g_training[:max_training_set_size]
			magnitude_r_training = magnitude_r_training[:max_training_set_size]
			magnitude_i_training = magnitude_i_training[:max_training_set_size]
			magnitude_z_training = magnitude_z_training[:max_training_set_size]
			magnitude_y_training = magnitude_y_training[:max_training_set_size]
			photo_z_training = photo_z_training[:max_training_set_size]
			galaxy_type_training = galaxy_type_training[:max_training_set_size]

		typecut_training = (galaxy_type_training < 9)

		print "... ... capped to max training set size"

		##   train

		# collect mags and photo-z in an array
#		M_training = numpy.array([magnitude_g_training,magnitude_r_training,magnitude_i_training,magnitude_z_training,magnitude_y_training,photo_z_training]).T
		color_gr_training = magnitude_g_training - magnitude_r_training
		color_ri_training = magnitude_r_training - magnitude_i_training
		color_iz_training = magnitude_i_training - magnitude_z_training
		color_zy_training = magnitude_z_training - magnitude_y_training
#		M_training = numpy.array([magnitude_g_training,magnitude_r_training,magnitude_i_training,magnitude_z_training,magnitude_y_training,color_gr_training,color_ri_training,color_iz_training,color_zy_training,photo_z_training]).T
		M_training = numpy.array([magnitude_g_training,magnitude_r_training,magnitude_i_training,magnitude_z_training,magnitude_y_training]).T
		numberOfColumns = M_training.shape[1]
		numberOfRows = M_training.shape[0]

		print "... ... collected data in M_training"

		# polynomial combinations
		for i in numpy.arange(0,numberOfColumns):
			for j in numpy.arange(i,numberOfColumns):
				M_training = numpy.hstack((M_training,(M_training[:,i]*M_training[:,j]).reshape(numberOfRows,1)))

		for i in numpy.arange(0,numberOfColumns):
			for j in numpy.arange(i,numberOfColumns):
				for k in numpy.arange(j,numberOfColumns):
					M_training = numpy.hstack((M_training,(M_training[:,i]*M_training[:,j]*M_training[:,k]).reshape(numberOfRows,1)))

#		for i in numpy.arange(0,numberOfColumns):
#			for j in numpy.arange(i,numberOfColumns):
#				for k in numpy.arange(j,numberOfColumns):
#					for l in numpy.arange(k,numberOfColumns):
#						M_training = numpy.hstack((M_training,(M_training[:,i]*M_training[:,j]*M_training[:,k]*M_training[:,l]).reshape(numberOfRows,1)))

		print "... ... added polynomial combinations, size: ", M_training.shape

		M_training = numpy.matrix(M_training)


		# Let's least-squares!
		b = numpy.matrix(typecut_training).T
		b = b.astype(float)
		[weights, residues, rank, singular_values] = numpy.linalg.lstsq(M_training,b)

#		pyplot.plot(b,M_training*weights,'.',b,b)
#		pyplot.show()

		b_estimated = (M_training*weights > .5)
		print "galaxies total:", magnitude_g_training.shape[0]
		print "should be selected:", numpy.sum(b)
		print "actually selected:", numpy.sum(b_estimated)
		print "correct positives:", (b.T*b_estimated)[0][0]
		print "false positives:", ((1-b).T*b_estimated)[0][0]
		print "correct negatives:", ((1 - b).T*(1 - b_estimated))[0][0]
		print "false negatives:", (b.T*(1-b_estimated))[0][0]


		print "... ... finished least squares, writing weights to disk"

		#   write weights to disk
		with archive.archive(weightsfile,'w') as ar:
			ar['/weights/weights'] = weights

		print "... finished training"



	# read weights from disk

	with archive.archive(weightsfile,'r') as ar:
		weights = ar['/weights/weights']


	print "... read weights from disk"

	# calculate cutoff function


	color_gr = magnitude_g - magnitude_r
	color_ri = magnitude_r - magnitude_i
	color_iz = magnitude_i - magnitude_z
	color_zy = magnitude_z - magnitude_y


	# collect mags and photo-z in an array
#	M = numpy.array([magnitude_g_cut,magnitude_r_cut,magnitude_i_cut,magnitude_z_cut,magnitude_y_cut,photo_z_cut]).T
#	M = numpy.array([magnitude_g,magnitude_r,magnitude_i,magnitude_z,magnitude_y,color_gr,color_ri,color_iz,color_zy,photo_z]).T
	M = numpy.array([magnitude_g,magnitude_r,magnitude_i,magnitude_z,magnitude_y]).T
	numberOfColumns = M.shape[1]
	numberOfRows = M.shape[0]

	print "... collected data in M"


	# polynomial combinations
	for i in numpy.arange(0,numberOfColumns):
		for j in numpy.arange(i,numberOfColumns):
			M = numpy.hstack((M,(M[:,i]*M[:,j]).reshape(numberOfRows,1)))

	for i in numpy.arange(0,numberOfColumns):
		for j in numpy.arange(i,numberOfColumns):
			for k in numpy.arange(j,numberOfColumns):
				M = numpy.hstack((M,(M[:,i]*M[:,j]*M[:,k]).reshape(numberOfRows,1)))

#	for i in numpy.arange(0,numberOfColumns):
#		for j in numpy.arange(i,numberOfColumns):
#			for k in numpy.arange(j,numberOfColumns):
#				for l in numpy.arange(k,numberOfColumns):
#					M = numpy.hstack((M,(M[:,i]*M[:,j]*M[:,k]*M[:,l]).reshape(numberOfRows,1)))


	M = numpy.matrix(M)

	print "... added polynomial combinations, size: ", M.shape

#	galaxy_type_estimated_cut = M*weights
#
#	galaxy_type_estimated = numpy.zeros(magnitude_g.shape)
#	galaxy_type_estimated[simplecuts_flag] = galaxy_type_estimated_cut

#	galaxy_type_estimated = numpy.array(M*weights)
	typecut_estimated = numpy.array(M*weights)


	print "... estimated galaxy type"

	# perform type cutoff

#	type_cutoff_flag = (galaxy_type_estimated < type_cutoff)*(galaxy_type_estimated > 0)
	# type_cutoff_flag = type_cutoff_flag[:,0] # convert from 2D (N-by-1) to 1D array
	type_cutoff_flag = (typecut_estimated > .5)

	print "... computed type-cutoff flag"



	# combine simple cuts and type cut
	target_selection_flag = numpy.empty(magnitude_g.shape,dtype=bool)
	target_selection_flag[:] = False
#	target_selection_flag[simplecuts_flag] = type_cutoff_flag
	target_selection_flag[simplecuts_flag] = type_cutoff_flag[simplecuts_flag]

	print "... computed target selection flag"

	# exporting selection flag
	with archive.archive(outputfile,'a') as ar:
		ar['/gal/target_selection_flag'] = target_selection_flag

	print "... wrote target selection flag to disk"




else:
	# use external training catalog to calculate the weights

	if (recalculate_weights):
		# import training cat, including selection flags
		training_catalog = csvimport.importArrayFromCSVFile(training_catalog_location)
		magnitude_g_training = training_catalog[:,0]
		magnitude_r_training = training_catalog[:,1]
		magnitude_i_training = training_catalog[:,2]
		magnitude_z_training = training_catalog[:,3]
		magnitude_y_training = training_catalog[:,4]
		photo_z_training = training_catalog[:,5]
		selection_flag_training = training_catalog[:,-1]

		# perform mag and photo-z cuts
		magcut_g_training_flag = (magnitude_g_training > mag_g_min)*(magnitude_g_training < mag_g_max)
		magcut_r_training_flag = (magnitude_r_training > mag_r_min)*(magnitude_r_training < mag_r_max)
		magcut_i_training_flag = (magnitude_i_training > mag_i_min)*(magnitude_i_training < mag_i_max)
		magcut_z_training_flag = (magnitude_z_training > mag_z_min)*(magnitude_z_training < mag_z_max)
		magcut_y_training_flag = (magnitude_y_training > mag_y_min)*(magnitude_y_training < mag_y_max)

		magcut_flag = magcut_g_training_flag*magcut_r_training_flag*magcut_i_training_flag*magcut_z_training_flag*magcut_y_training_flag

		photozcut_training_flag = (photo_z_training > photo_z_min)*(photo_z_training < photo_z_max)

		simplecuts_training_flag = magcut_training_flag*photozcut_training_flag

		magnitude_g_training = magnitude_g_training[simplecuts_training_flag]
		magnitude_r_training = magnitude_r_training[simplecuts_training_flag]
		magnitude_i_training = magnitude_i_training[simplecuts_training_flag]
		magnitude_z_training = magnitude_z_training[simplecuts_training_flag]
		magnitude_y_training = magnitude_y_training[simplecuts_training_flag]
		photo_z_training = photo_z_training[simplecuts_training_flag]
		galaxy_type_training = galaxy_type_training[simplecuts_training_flag]

		print "... ...  performed simple cuts on training cat"

		#   cap to max_training_set_size
		if (numpy.sum(simplecuts_training_flag) > max_training_set_size):
			magnitude_g_training = magnitude_g_training[:max_training_set_size]
			magnitude_r_training = magnitude_r_training[:max_training_set_size]
			magnitude_i_training = magnitude_i_training[:max_training_set_size]
			magnitude_z_training = magnitude_z_training[:max_training_set_size]
			magnitude_y_training = magnitude_y_training[:max_training_set_size]
			photo_z_training = photo_z_training[:max_training_set_size]
			galaxy_type_training = galaxy_type_training[:max_training_set_size]

		print "... ... capped to max training set size"



		# collect mags and photo-z in an array
		M_training = numpy.array([magnitude_g_training,magnitude_r_training,magnitude_i_training,magnitude_z_training,magnitude_y_training,photo_z_training]).T
		numberOfColumns = M_training.shape[1]
		numberOfRows = M_training.shape[0]

		print "... ... collected data in M_training"

		# polynomial combinations
		for i in numpy.arange(0,numberOfColumns):
			for j in numpy.arange(i,numberOfColumns):
				M_training = numpy.hstack((M_training,(M_training[:,i]*M_training[:,j]).reshape(numberOfRows,1)))

		print "... ... added polynomial combinations, size: ", M_training.shape

		M_training = numpy.matrix(M_training)


		# Let's least-squares!
		b = numpy.matrix(selection_flag_training).T
		b = b.astype(float)
		[weights, residues, rank, singular_values] = numpy.linalg.lstsq(M_training,b)


		print "... ... finished least squares, writing weights to disk"

		#   write weights to disk
		with archive.archive(inputfile,'a') as ar:
			ar['/target_selection_training/weights_external'] = weights

		print "... finished training"


	# read weights from disk

	with archive.archive(outputfile,'r') as ar:
		weights = ar['/target_selection_training/weights_external']

	print "... read weights from disk"

	# calculate cutoff function

	# collect mags and photo-z in an array
	M = numpy.array([magnitude_g,magnitude_r,magnitude_i,magnitude_z,magnitude_y,photo_z]).T
	numberOfColumns = M.shape[1]
	numberOfRows = M.shape[0]

	print "... collected data in M"


	# polynomial combinations
	for i in numpy.arange(0,numberOfColumns):
		for j in numpy.arange(i,numberOfColumns):
			M = numpy.hstack((M,(M[:,i]*M[:,j]).reshape(numberOfRows,1)))

	M = numpy.matrix(M)

	print "... added polynomial combinations, size: ", M.shape

	selection_flag_estimated_float = M*weights

	print "... estimated selection flag"

	# perform cutoff
	# for that we need to estimate where to cut off between 0 and 1
	below_cutoff = numpy.sum((selection_flag_training == 0))
	above_cutoff = numpy.sum((selection_flag_training == 1))
	cutoff = below_cutoff/(below_cutoff + above_cutoff)

	selection_flag_estimated = (selection_flag_estimated_float > cutoff)

	target_selection_flag = selection_flag_estimated




print "... computed target selection flag"

numberOfGalaxies = magnitude_g.size
selected = numpy.sum((target_selection_flag == True))
discarded = numpy.sum((target_selection_flag == False))
undefined = magnitude_g.size - selected - discarded
print "... selected:", selected
print "... discarded:", discarded
print "... undefined:", undefined
print "=== total:", numberOfGalaxies
print
print "****************************"
if (undefined == 0):
	print "*** checksum passed ***"
else:
	print "XXX CHECKSUM FAILED XXX"
print "****************************"


# exporting selection flag
with archive.archive(outputfile,'a') as ar:
	ar['/gal/target_selection_flag'] = target_selection_flag

print "... wrote target selection flag to disk"


# Diagnostic
# gmr vs. rmi
# gmr vs. rmag
# gmr vs. z_true
print "... plotting diagnostics now ..."

numberOfGalaxiesPlotted = 1000

discarded_flag = numpy.logical_not(target_selection_flag)*magcut99_flag

numberOfGalaxiesWithout99s = numpy.sum(magcut99_flag)
numberOfGalaxiesSelected = numpy.sum(target_selection_flag)
numberOfGalaxiesDiscarded = numberOfGalaxiesWithout99s - numberOfGalaxiesSelected

gmag_sel = magnitude_g[target_selection_flag]
gmag_dis = magnitude_g[discarded_flag]
rmag_sel = magnitude_r[target_selection_flag]
rmag_dis = magnitude_r[discarded_flag]
imag_sel = magnitude_i[target_selection_flag]
imag_dis = magnitude_i[discarded_flag]
zmag_sel = magnitude_z[target_selection_flag]
zmag_dis = magnitude_z[discarded_flag]
ymag_sel = magnitude_y[target_selection_flag]
ymag_dis = magnitude_y[discarded_flag]
photoz_sel = photo_z[target_selection_flag]
photoz_dis = photo_z[discarded_flag]
#type_selected = galaxy_type_estimated[target_selection_flag]
#type_discarded = galaxy_type_estimated[discarded_flag]



#plot_type='scatter'
plot_type='hist'
show_plot=False

# g-r vs r
Plotters.plot_target(rmag_sel,gmag_sel-rmag_sel ,x2=rmag_dis, y2=gmag_dis-rmag_dis,
					xmin=19, xmax=29,
					ymin=-1, ymax=3.5,
					plot_type=plot_type,
					Nbinx=10, Nbiny=10,
					#cmap_name1=cmap_name1,
					label_x='r mag', label_y='(g-r)',
					fig_number=0, plot_number=0, title='color-mag of targets', base_directory='.', filename_addendum="gmr_v_r"+"_",
					show_plot=show_plot)


# r-i vs i
Plotters.plot_target(imag_sel,rmag_sel-imag_sel ,x2=imag_dis, y2=rmag_dis-imag_dis,
					xmin=19, xmax=29,
					ymin=-1, ymax=3.5,
					plot_type=plot_type,
					Nbinx=10, Nbiny=10,
					label_x='i mag', label_y='(r-i)',
					fig_number=1, plot_number=0, title='color-mag of targets', base_directory='.', filename_addendum="gmr_v_g"+"_",
					show_plot=show_plot)

# g-r vs r-i
Plotters.plot_target(rmag_sel-imag_sel,gmag_sel-rmag_sel ,x2=rmag_dis-imag_dis, y2=gmag_dis-rmag_dis,
					xmin=-1, xmax=3.5,
					ymin=-1, ymax=3.5,
					plot_type=plot_type,
					Nbinx=10, Nbiny=10,
					label_x='(r-i)', label_y='(g-r)',
					fig_number=2, plot_number=0, title='color-color of targets', base_directory='.', filename_addendum="gmr_v_rmi"+"_",
					show_plot=show_plot)

# g-r vs z
Plotters.plot_target(photoz_sel, gmag_sel-rmag_sel ,x2=photoz_dis, y2=gmag_dis-rmag_dis,
					xmin=0, xmax=2.,
					ymin=-1, ymax=3.5,
					plot_type=plot_type,
					Nbinx=10, Nbiny=10,
					label_x=r'$z_{photo}$', label_y='(g-r)',
					fig_number=3, plot_number=0, title='color-z of targets', base_directory='.', filename_addendum="gmr_v_z"+"_",
					show_plot=show_plot)

# r-i vs z
Plotters.plot_target(photoz_sel, rmag_sel-imag_sel ,x2=photoz_dis, y2=rmag_dis-imag_dis,
					xmin=0, xmax=2.,
					ymin=-1, ymax=3.5,
					plot_type=plot_type,
					Nbinx=10, Nbiny=10,
					label_x=r'$z_{photo}$', label_y='(r-i)',
					fig_number=4, plot_number=0, title='color-zof targets', base_directory='.', filename_addendum="rmi_v_z"+"_",
					show_plot=show_plot)


#pyplot.xlabel('i')
#pyplot.ylabel('type')
#pyplot.ylim((-5,20))
#pyplot.title('Target Selection: mag-type')
#p1 = pyplot.scatter(imag_dis,type_discarded,marker='o',color='r')
#p2 = pyplot.scatter(imag_sel,type_selected,marker='D',color='b')
#pyplot.legend([p1,p2],["discarded ("+str(numberOfGalaxiesDiscarded)+")","selected ("+str(numberOfGalaxiesSelected)+")"])
#pyplot.savefig(file_figure_out_i_type)










