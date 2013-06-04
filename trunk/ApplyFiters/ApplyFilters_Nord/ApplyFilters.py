# import
import sys

# import user-def functions
sys.path.append('../../Wrapper/')
import archive
sys.path.append('../../Plotters/')
import Plotters

# import
import numpy



# ================================================
# Detect Quantity
# ================================================
def GenerateFilters(wavelengths, wavelength_fiducial, resolution_fiducial, wavelength_min, wavelength_max):

    # FWHM model
    fwhm = lambda x: x**2/ (wavelength_fiducial* resolution_fiducial)

    # calculate deviatison from quantities
    fwhms = numpy.array([ numpy.array(fwhm(lambdaw)) for lambdaw in wavelengths])

    return fwhms


# ================================================
# Generate Photons
# ================================================
def DetectPhotons():

    # -----------------------------------------------------
    # Import data and params
    # -----------------------------------------------------
    print
    print "........................."
    print " Entering Apply Filters "
    print "........................."

    databank_file = '../../../data/data_bank.h5'
    with archive.archive(databank_file, 'r') as ar:
        # filter information
        wavelength_fiducial = ar['/Instrument/Filter/wavelength_fiducial']
        resolution_fiducial = ar['/Instrument/Filter/resolution_fiducial']
        wavelength_min      = ar['/Instrument/Filter/wavelength_minimum']
        wavelength_max      = ar['/Instrument/Filter/wavelength_maximum']

        # photons
        spectrum            = ar[]

    # -----------------------------------------------------
    # Generate Filters
    # -----------------------------------------------------
    Filters = GenerateFilters(wavelengths, wavelength_fiducial, resolution_fiducial)

    # -----------------------------------------------------
    # Apply Filters
    # ... later: apply filter to each spectrum
    # -----------------------------------------------------
# loop over galaxies
    for gal in
    # ... apply filter to photon list how are photons used in EAZY photo-z calculator


    # -----------------------------------------------------
    # Export data and params
    # -----------------------------------------------------
    with archive.archive(databank_file, 'r') as ar:

    print
    print "........................."
    print " Exiting Generate Image  "
    print "........................."


    # plot diagnostics
    # do the following for all new photon information
    #   overlay photon energies on initial galaxy spectrum
    #   plot photon times
    #   positions vs. time
    #   positions vs. energy
    #   energy vs. time
    # compare all old vs. new photon information


    # Quality Assurance




# ================================================
# Main
# ================================================
def main():

    ApplyDetectorMeasurement()




#-------------------------------
if __name__ == '__main__':
    main()
