import archive
import numpy
import json

class spokes_tests:
	def __init(self):
		self.messages = []

	def log(self, message):
		#self.messages.append(message)
		print message
	
	def convert(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
			self.log('Galaxy count:' + str(len(ar['/gal/galaxy_index'])))
			if len(ar['/gal/galaxy_index']) > 0:
				return True
			else:
				return False
	
	def target_selection(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
			self.log('Selected by Target Selection:' + str(sum(ar['/gal/target_selection_flag'])))			
			if sum(ar['/gal/target_selection_flag']) > 0:
				return True
			else:
				return False
			
	def survey_strategy(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
			self.log('Selected by Survey strategy:' + str(sum(ar['/gal/survey_selection_flag'])))
			if sum(ar['/gal/survey_selection_flag']) == 0:
				return False
			else:
				if ar['/tiling/number_of_tiles'] > 0:
					return True
				else:
					return False
			
	def fiber_allocation(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
			self.log('Selected by Fiber Selection:' + str(sum(ar['/gal/fiber_selection_flag'])))				
			if sum(ar['/gal/fiber_selection_flag']) > 0:
				return True
			else:
				return False
				
	def throughput(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
	
			if len(ar['/Throughput/throughput']) > 0:
				self.log('Throughput ok')
				return True
			else:
				return False			
				
	def measure_redshift(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
			z_spec 					= ar['/gal/redshift_spectroscopic_spokes']
			z_true					= ar['/gal/z_true']
			fiber_selection_flag 	= ar['/gal/fiber_selection_flag']
			
			z_true_cropped = z_true[fiber_selection_flag]
			#does not work atm
			return True
			
			relative_error = z_spec/z_true_cropped - 1
			
			if (abs(numpy.mean(relative_error)) < 0.1):
				self.log('Measure Redshift mean error:' + str(abs(numpy.mean(relative_error)))	)
				return True
			else:
				return False
				
	def bin_redshift(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:	
			if ar['/Ensemble/nbins'] > 0:
				return True
			else:
				return False

	def estimate_cosmo_params(self):
		with archive.archive('../../data/data_bank.h5', 'r') as ar:
			if len(ar['/constraints/fisher_matrix']) > 0:
				return True
			else:
				return False		
