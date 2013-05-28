# import
import sys

sys.path.append('../../Wrapper/')
import archive
sys.path.append('../../Plotters/')
import Plotters

import numpy
from scipy import stats


# ================================================
# Count Number of Photons
# ================================================
def magphoton(mag,exptime=1.,gain=0.23):            # !!!PAR
    """
        This codes generate the number of photons and ADUs based on DES optics and detectors.
    """
    zeropoint = 26.794176                           # !!!PAR
    nphoton = exptime*10**(0.4*(zeropoint - mag))
    nadu = nphoton*gain

    return nphoton,nadu



# ================================================
# Sample Position Domain
# ================================================
def SamplePositions(image, nphoton, seed_position=9999 ):  # !!!PAR

    # set random seed
    random.seed(seed_position)

    # reshape image into 1D array
        # track pixel locations
    # remove zero-value pixels
    # sample photon positions along the array (for each photon)
    # get positions from which pixels are sampled

    return pos_x, pos_y



# ================================================
# Sample Time Domain
# ================================================
def SampleTimeArrival(nphoton, exptime, start_time=0.0, seed_time=9999 ):  # !!!PAR

    # set random seed
    random.seed(seed_time)

    # sample from exposure time
    photon_arrival_times = [ random.uniform(start_time, exptime) for i in nphoton ]

    return photon_arrival_times



# ================================================
# Sample Energy Domain
# ================================================
def SampleEnergySpectrum(nphoton, seed_energy=9999 ):               # !!!PAR

    # set random seed
    random.seed(seed_energy)

    # set up distribution from galaxy spectrum
    spectrum_x = np.arange(len(spectrum))
    spectrum_sampler= stats.rv_discrete(name='spectrum', values=(spectrum_x, spectrum))
    #h = plt.plot(spectrum_x, spectrum_sampler.pmf(spectrum_x)) # check the result with this plot

    # sample from spectrum
    photon_energies = spectrum_sampler.rvs(size=nphoton)

    return photon_energies



# ================================================
# Generate Photons
# ================================================
def GeneratePhotons():

    # -----------------------------------------------------
    # Import data and params
    # -----------------------------------------------------
    print
    print "........................."
    print " Entering Generate Image  "
    print "........................."

    databank_file = '../../../data/data_bank.h5'
    with archive.archive(databank_file, 'r') as ar:
        exptime              = ar['/SurveyDesign/exposure_time']
        magnitude_galaxy     = ar['/Galaxy/magnitude/']
        flux_spectrum_galaxy = ar['/Galaxy/spectrum']
        image                = ar['/Galaxy/image']


    # Calculate number of photons
    nphoton, = magphoton(magnitude_galaxy, exptime=exptime)


    # Sample Photon Positions
    PhotonPositionsX, PhotonPositionsY = SamplePositions(image, nphoton)


    # Sample Photon Times
    PhotonArrivalTimes = SampleTimeArrival(nphoton, exptime, start_time=0.0, seed_time=9999 )


    # Sample Photon Energies
    PhotonEnergies = SampleEnergySpectrum(nphoton)


    # Generate Photon Identification Numbers
    photon_id = numpy.arange(nphoton)


    # -----------------------------------------------------
    # Export data and params
    # -----------------------------------------------------
    with archive.archive(databank_file, 'r') as ar:
        ar['/Photon/Energy']            = PhotonEnergies
        ar['/Photon/ArrivalTime']       = PhotonArrivalTimes
        ar['/Photon/PositionsX']        = PhotonPositionsX
        ar['/Photon/PositionsY']        = PhotonPositionsY
        ar['/Photon/Number_of_Photons'] = nphoton

    print
    print "........................."
    print " Exiting Generate Image  "
    print "........................."


    # plot diagnostics
    #   overlay photon energies on initial galaxy spectrum
    #   plot photon times
    #   positions vs. time
    #   positions vs. energy
    #   energy vs. time

    # Quality Assurance



# ================================================
# Main
# ================================================
def main():

    GeneratePhotons()




#-------------------------------
if __name__ == '__main__':
    main()
