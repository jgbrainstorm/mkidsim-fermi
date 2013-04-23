#---this file define the functions used in the simulation ---
# created by J. Hao, April 16, 2013
from scipy.misc import factorial as fac


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


def zernike_rad(m, n, rho):
    """
    Calculate the radial component of Zernike polynomial (m, n) 
    given a grid of radial coordinates rho.
    """
    if (n < 0 or m < 0 or abs(m) > n):
        raise ValueError
    if ((n-m) % 2):
        return rho*0.0
    pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
    return sum(pre_fac(k) * rho**(n-2.0*k) for k in xrange((n-m)/2+1))

def zernike(m, n, rho, phi):
    """
    Calculate Zernike polynomial (m, n) given a grid of radial
    coordinates rho and azimuthal coordinates phi.
    """
    if (m > 0): return zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0): return zernike_rad(-m, n, rho) * np.sin(-m * phi)
    return zernike_rad(0, n, rho)

def zernikel(j, rho, phi):
    """
    Calculate Zernike polynomial with Noll coordinate j given a grid of radial coordinates rho and azimuthal coordinates phi.
    """
    n = 0
    while (j > n):
        n += 1
        j -= n
    m = -n+2*j
    return zernike(m, n, rho, phi)


def dispZernike(beta=1.,j=0,gridsize = 1, max_rad = 1):
    x,y = np.meshgrid(np.arange(-gridsize,gridsize,0.001),np.arange(-gridsize,gridsize,0.001))
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    ok = rho < max_rad
    znk = beta*zernikel(j,rho,phi)*ok
    pl.imshow(znk)
    return znk
    

