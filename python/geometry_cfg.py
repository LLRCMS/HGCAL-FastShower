
import math as m

geometry_type = 'Hexagons'

#rotation angle of cells between two layers (degree)
#implemented for hexagons cells only
geometry_cell_rotation = 0.

geometry_layer = -1
# layers' z positions
# from CMSSW V7 geometry: https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
# the values are the silicon (centre) z positions of the 40 layers wrt HGCAL z entry position in cm
geometry_layers_z = [.75, 1.5, 2.73, 3.48, 4.71, 5.46, 6.69, 7.44, 8.67, 9.42, 10.73, 11.6, 12.91,
                     13.78, 15.09, 15.96, 17.27, 18.14, 19.45, 20.32, 21.77, 22.84, 24.29, 25.36,
                     26.81, 27.88, 29.33, 30.4, 36.33, 41.01, 45.69, 50.37, 55.05, 59.73, 64.41,
                     69.09, 73.77, 78.45, 83.13, 87.81, 105.1, 113.0, 122.5, 131.2, 139.9, 148.6,
                     157.3, 166.0, 174.7, 183.4, 192.1, 200.8]
# entry point in HGCal
z0 = 320.
# absolute HGCal layer positions
geometry_layers_z[:] = [z+z0 for z in geometry_layers_z]

# multiply cell side by sqrt(6) for triangles to get equal area
geometry_small_cell_side = 0.476
geometry_large_cell_side = 0.648

# limits position of the different zones (100, 200, 300 um)
if geometry_layer in range(0, 27): #EE
    geometry_limit_first_zone = 75. #cm
    geometry_limit_second_zone = 120.
else: #FH-BH
    geometry_limit_first_zone = 60. #cm
    geometry_limit_second_zone = 100.


# Define (eta,phi) window to build the geometry
geometry_eta_min = 1.475
geometry_eta_max = 3.
geometry_phi_min = -m.pi/8.
geometry_phi_max = m.pi/8.


geometry_file = 'data/AC_11.json'
