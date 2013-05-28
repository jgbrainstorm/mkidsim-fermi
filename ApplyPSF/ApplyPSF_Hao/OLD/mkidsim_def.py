#---this file define the functions used in the simulation ---
# created by J. Hao, April 16, 2013


from math import *
import numpy as np
import numpy.fft
import sys,time
from scipy.special import gamma
from numpy.random import normal
import scipy.fftpack as fft



def magphoton(mag,exptime=1.,gain=0.23):
    """
    This codes generate the number of photons and ADUs based on DES optics and detectors. 
    """
    zeropoint = 26.794176
    nphoton = exptime*10**(0.4*(zeropoint - mag))
    nadu = nphoton*gain
    return nphoton,nadu



def des_psf_image(exptime=100,mag=None,seeing=[0.7,0.,0.],setbkg=True,moffat=False):
    
    """
    This code generate a PSF star with seeing and sky background (no optics psf)
    exptime is given in sec
    seeing is give in terms of [fwhm (arcsec),e1,e2]
    """
    gain = 0.21 # convert electrons to ADU
    npix = 40
    zeropoint = 26.794176 # r band, from Nikolay
    objectphoton = exptime*10**(0.4*(zeropoint - mag))
    if setbkg == False:
        skyphoton = 0.
    else:
        skyphoton = 8.460140*exptime #(sky level per pix per sec)
    bkg = skyphoton*gain  # background in ADU
    if moffat == True:
        psf = moffat_psf(npix = npix,fwhm=seeing[0],beta=3.5,scale=0.27)
    else:
        psf = gauss_seeing(npix,seeing[0],seeing[1],seeing[2],scale = 0.27)
    img = (psf * objectphoton + skyphoton)*gain
    img = img + add_imageNoise(img)
    return img,bkg,psf


def plane_wavefront(npix=256):
    wf = np.zeros((npix,npix),dtype='d')
    return wf



def factorial(n):
    return gamma(n+1)



def zernike(j,npix=256,phase=0.0):
    '''
    generate zernike polynomial of jth order
    '''
    if (j > 820):
      print "For n < 40, pick j < 820"
      sys.exit()
    x = np.arange(-npix/2,npix/2,dtype='d')
    y = np.arange(-npix/2,npix/2,dtype='d')
    xarr = np.outer(np.ones(npix,dtype='d'),x)
    yarr = np.outer(y,np.ones(npix,dtype='d'))
    rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
    thetarr = np.arctan2(yarr,xarr) + phase
    outside = np.where(rarr > 1.0)
    narr = np.arange(40)
    jmax = (narr+1)*(narr+2)/2
    wh = np.where(j <= jmax)
    n = wh[0][0]
    mprime = j - n*(n+1)/2
    if ((n % 2) == 0):
      m = 2*int(floor(mprime/2))
    else:
      m = 1 + 2*int(floor((mprime-1)/2))
    radial = np.zeros((npix,npix),dtype='d')
    zarr = np.zeros((npix,npix),dtype='d')
    for s in range((n-m)/2 + 1):
      tmp = pow(-1,s) * factorial(n-s)
      tmp /= factorial(s)*factorial((n+m)/2 - s)*factorial((n-m)/2 - s)
      radial += tmp*np.power(rarr,n-2*s)
    if (m == 0):
      zarr = radial
    else:
      if ((j % 2) == 0):
        zarr = sqrt(2.0)*radial*np.cos(m*thetarr)
      else:
        zarr = sqrt(2.0)*radial*np.sin(m*thetarr)
    zarr *= sqrt(n+1)
    zarr[outside] = 0.0
    return zarr


def aperture(npix=256, cent_obs=0.0, spider=0):
    '''
    this function make a spide like mask
    '''
    illum = np.ones((npix,npix),dtype='d')
    x = np.arange(-npix/2,npix/2,dtype='d')
    y = np.arange(-npix/2,npix/2,dtype='d')
    xarr = np.outer(np.ones(npix,dtype='d'),x)
    yarr = np.outer(y,np.ones(npix,dtype='d'))
    rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
    outside = np.where(rarr > 1.0)
    inside = np.where(rarr < cent_obs)
    illum[outside] = 0.0
    if np.any(inside[0]):
        illum[inside] = 0.0
    if (spider > 0):
        start = npix/2 - int(spider)/2
        illum[start:start+int(spider),:] = 0.0
        illum[:,start:start+int(spider)] = 0.0
    return illum


def wavefront(d_over_r0, npix=256, nterms=15, level=None):
    '''
    this function generate the phase (wavefront) by assuming Kolmogorov turbulence model.
    a description of this model can be fround from: 
    http://www.ctio.noao.edu/~atokovin/tutorial/part1/turb.html
    '''
    scale = pow(d_over_r0,5.0/3.0)
    if (nterms < 10):
        print "use at least use ten terms..."
        sys.exit()
    if level:
        narr = np.arange(400,dtype='d') + 2
        coef = np.sqrt(0.2944*scale*(np.power((narr-1),-0.866) - np.power(narr,-0.866)))
        wh = np.where(coef < level)
        n = wh[0][0]
        norder = int(ceil(sqrt(2*n)-0.5))
        nterms = norder*(norder+1)/2
        if (nterms < 15):
            nterms = 15
    wf = np.zeros((npix,npix),dtype='d')
    resid = np.zeros(nterms,dtype='d')
    coeff = np.zeros(nterms,dtype='d')
    resid[0:10] = [1.030,0.582,0.134,0.111,0.088,0.065,0.059,0.053,0.046,0.040]
    if (nterms > 10):
        for i in range(10,nterms):
            resid[i] = 0.2944*pow(i+1,-0.866)
    for j in range(2,nterms+1):
        coeff[j-1] = sqrt((resid[j-2]-resid[j-1])*scale)
        wf += coeff[j-1]*normal()*zernike(j,npix=npix)
    return wf


def atmPSF(aperture, wavefront, overfill=1):
    '''
    this codes generate the PSF due to atmospheric turbulence. 
    '''
    npix = len(wavefront)
    nbig = npix*overfill
    wfbig = np.zeros((nbig,nbig),dtype='d')
    half = (nbig - npix)/2
    wfbig[half:half+npix,half:half+npix] = wavefront
    illum = np.zeros((nbig,nbig),dtype='d')
    illum[half:half+npix,half:half+npix] = aperture
    phase = np.exp(wfbig*(0.+1.j))
    input = illum*phase
    ft = fft.fft2(input)
    ft = fft.fftshift(ft)
    powft = abs(ft)**2
    crop =  powft[half:half+npix,half:half+npix]
    return crop


def genPSF(d_over_r0=20.,nterms=20,overfill=1.5):
    aper = aperture(npix=256, cent_obs=0.2, spider=2)
    wf = wavefront(d_over_r0, npix=256, nterms=nterms, level=None)
    img = atmPSF(aper,wf,overfill)
    #pl.imshow(img)
    return img/img.sum()


    
for i in range(100):
    im = genPSF()
    pl.imshow(im)
    pl.savefig(str(i)+'.png')
    pl.close()
    print i



