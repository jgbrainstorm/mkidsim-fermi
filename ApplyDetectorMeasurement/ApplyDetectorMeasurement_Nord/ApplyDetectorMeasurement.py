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
def DetectQuantity(quantity, sigma, seed_quantity=9999):

    # set random seed
    random.seed(seed_quantity)

    # calculate deviatison from quantities
    quantity_delta = [ random.normalvariate(quantity[i], sigma) for i in len(quantity)]

    # add deviations to quantities
    quantity += quantity_delta

    return quantity


# ================================================
# Generate Photons
# ================================================
def DetectPhotons():

    # -----------------------------------------------------
    # Import data and params
    # -----------------------------------------------------
    print
    print "........................."
    print " Entering Generate Image  "
    print "........................."

    databank_file = '../../../data/data_bank.h5'
    with archive.archive(databank_file, 'r') as ar:
        PhotonEnergies      = ar['/Photon/Energy']
        PhotonArrivalTimes  = ar['/Photon/ArrivalTime']
        PhotonPositionsX    = ar['/Photon/PositionsX']
        PhotonPositionsY    = ar['/Photon/PositionsY']

        position_resolution = ar['/Instrument/Detector/position_resolution']
        time_resolution     = ar['/Instrument/Detector/time_resolution']
        energy_resolution   = ar['/Instrument/Detector/energy_resolution']

        seed_position_x     = 9999
        seed_position_y     = 9999
        seed_time           = 9999
        seed_energy         = 9999


    # Detect Photon Positions
    PhotonPositionsX_Detect = DetectQuantity(PhotonPositionsX, position_resolution, seed_quantity=seed_position_x)
    PhotonPositionsY_Detect = DetectQuantity(PhotonPositionsY, position_resolution, seed_quantity=seed_position_y)


    # Detect Photon Times
    PhotonArrivalTimes_Detect = DetectQuantity(PhotonArrivalTimes, time_resolution, seed_quantity=seed_time)


    # Detect Photon Energies
    PhotonEnergies_Detect = DetectQuantity(PhotonEnergies, position_resolution, seed_quantity=seed_energy)


    # -----------------------------------------------------
    # Export data and params
    # -----------------------------------------------------
    with archive.archive(databank_file, 'r') as ar:
        ar['/Photon/Energy_Detected']       = PhotonEnergies
        ar['/Photon/ArrivalTime_Detected']  = PhotonArrivalTimes
        ar['/Photon/PositionsX_Detected']   = PhotonPositionsX
        ar['/Photon/PositionsY_Detected']   = PhotonPositionsY
        ar['/Photon/Number_of_Photons']     = nphoton

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
