
import math as m

#choose E gun or pT gun
generation_gun_type = "E"

# Energy is a generic term, enter your E or pT values
generation_energy = 35. #GeV
generation_energy_fluctuation = False
generation_energy_range_max = 120.
generation_energy_range_min = 100.

# Incident direction
generation_incident_eta = 2.
generation_eta_fluctuation = False
generation_eta_range_max = 2.9
generation_eta_range_min = 3.

generation_incident_phi = 0.
generation_phi_fluctuation = False
generation_phi_range_max = 0.01
generation_phi_range_min = -0.01


generation_fluctuation = True
# nbr hits per GeV, used only if fluctuation is false
generation_number_of_hits_per_gev = 1000


# electronic noise
generation_noise = True
# Calibration (noise true) : you can choose to calibrate with values from file or with your own values
generation_calib_type = "External"
generation_file = './calibration.txt'

# mev per MIP for 200 um
generation_mip = [9.4, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 9.7, 11, 11, 11, 11, 11, 11,
                       11, 11, 11, 12.7, 14.8, 14.8, 14.8, 14.8, 14.8, 14.8, 14.8, 41.5, 59.8, 54.2,
                       54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 67.4, 83.2, 92.2, 92.2,
                       92.2, 92.2, 92.2, 92.2, 92.2, 92.2, 92.2, 92.2, 46.1]

 # noise in mips for 200 um
generation_noise_sigma = [0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
