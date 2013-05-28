# import
import sys
sys.path.append('../../Wrapper/')
sys.path.append('../../Plotters/')
import archive
import Plotters
import numpy


# ================================================
# Generate Multivariate Normal
# ================================================
def GenerateMultivariateNormal(pos_x, pos_y, center_x, center_y, sigma_x, sigma_y, correlation):
    A = 1./ (2.* numpy.pi* sigma_x* sigma_y* numpy.sqrt(1-corr**2) )
    B0 = -1./ (2.* (1-corr**2))
    B1 = (pos_x - center_x)**2/ sigma_x**2
    B2 = (pos_y - center_y)**2/ sigma_y**2
    B3 = 2.* corr* (pos_x - center_x)* (pos_y - center_y)/ (sigma_x* sigma_y)
    B = B0* (B1 + B2 + B3)
    return A* B



# ================================================
# Measure Center
# ================================================
def MeasureCenter(npix):
    return numpy.ceil(float(npix)/2.-1.), numpy.ceil(float(npix)/2.-1.)



# ================================================
# Generate Profile
# ================================================
def GenerateProfile(image, npix, profiletype):
    center_x, center_y = MeasureCenter(npix)
    if profiletype == 'center':
        image[center_x,center_y] = 1.
    if profiletype == 'gaussian':
        print "does not exist yet"



# ================================================
# Generate Image
# ================================================
def GenerateImage():

    # -----------------------------------------------------
    # Import data and params
    # -----------------------------------------------------
    print
    print "........................."
    print " Entering Generate Image  "
    print "........................."

    databank_file = '../../../data/data_bank.h5'
    with archive.archive(databank_file, 'r') as ar:
        npix = ar['/Instrument/Detector/Number_of_Pixels']
        profile_type = ar['/Instrument/Detector/Profile_Type']


    # -----------------------------------------------------
    # Perform Calculations
    # -----------------------------------------------------

    # Create 2-D Array of points
    image = numpy.zeros( (npix,npix), dtype=float)

    # create profile and add to image
    image = GenerateProfile( image, npix, profile_type )


    # -----------------------------------------------------
    # Export data and params
    # -----------------------------------------------------
    with archive.archive(databank_file, 'r') as ar:
        ar['/Galaxy/image']

    print
    print "........................."
    print " Exiting Generate Image  "
    print "........................."


    # -----------------------------------------------------
    # Plotting Diagnostics
    # -----------------------------------------------------
# plot diagnostics
#   image
#   1d profile plot
#   histogram of pixel brightness


    # -----------------------------------------------------
    # Quality Assessment
    # -----------------------------------------------------





# ================================================
# Main
# ================================================
def main():

    GenerateImage()




#-------------------------------
if __name__ == '__main__':
    main()
