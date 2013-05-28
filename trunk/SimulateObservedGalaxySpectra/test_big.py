# Copyright 2012 ETHZ.ch Lukas Gamper, gamperl@gmail.com

import numpy as np
from kcorrect import k_binspec
import random

nbin_old = 100000
nbin_new = 200000

lambdaArr   = np.array( range(nbin_old),  dtype='float32')
newLambdaArr   = np.array( range(nbin_new),  dtype='float32')
spectrumArr = np.array( [random.randint(1,1000) for x in xrange(nbin_old)],  dtype='float32')

newSpectrumArr = np.array(np.zeros(newLambdaArr.shape), dtype='float32')

k_binspec(lambdaArr, spectrumArr, newLambdaArr, newSpectrumArr)

print newSpectrumArr
