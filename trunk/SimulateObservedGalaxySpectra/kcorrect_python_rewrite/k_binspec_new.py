def k_binspec(lambda_arr, spectrum, newlambda_arr, newspectrum):
	lambda_diff = lambda_arr[1:] - lambda_arr[:-1]
	spectrum_diff = (spectrum[:-1] + spectrum[1:]) * .5
	i = 0
	j = 0
	while i < len(lambda_arr) - 1 and j < len(newlambda_arr):
		while lambda_arr[i + 1] < newlambda_arr[j]:
			newspectrum[j] += lambda_diff[i] * spectrum_diff[i]
			i += 1
		while i < len(lambda_arr) - 1 and j < len(newlambda_arr) and lambda_arr[i + 1] > newlambda_arr[j]:
			newspectrum[j] += (newlambda_arr[j] - lambda_arr[i]) * spectrum_diff[i]
			newspectrum[j] += (lambda_arr[i + 1] - newlambda_arr[j]) * spectrum_diff[i]
			j += 1
		i += 1
