# Copyright 2012 ETHZ.ch Lukas Gamper, gamperl@gmail.com

import numpy as np
from kcorrect import k_binspec


lambdaArr = np.array([1, 2, 3, 4, 5], dtype='float32')
spectrumArr = np.array([1, 2, 3, 4, 5], dtype='float32')
newLambdaArr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype='float32')
newSpectrumArr = np.array(np.zeros(newLambdaArr.shape), dtype='float32')

k_binspec(lambdaArr, spectrumArr, newLambdaArr, newSpectrumArr)

print newSpectrumArr
