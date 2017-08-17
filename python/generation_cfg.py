
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
generation_calib_type = "Personnal"
generation_file = 'calib.txt'

# # mev per MIP (GeV) for 200 um
generation_mip_energy = [0.094, 0.088, 0.088, 0.088, 0.088, 0.088, 0.088, 0.088, 0.088, 0.097,
                        0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.127, 0.148, 
                        0.148, 0.148, 0.148, 0.148, 0.148, 0.148, 0.415, 0.598, 0.542, 0.542,
                        0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.674, 0.832, 0.922,
                        0.922, 0.922, 0.922, 0.922, 0.922, 0.922, 0.922, 0.922, 0.922, 0.461]
# sampling fraction for 200 um
generation_sampling = [9.4, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 9.7, 11, 11, 11, 11, 11, 11,
                       11, 11, 11, 12.7, 14.8, 14.8, 14.8, 14.8, 14.8, 14.8, 14.8, 41.5, 59.8, 54.2,
                       54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 54.2, 67.4, 83.2, 92.2, 92.2,
                       92.2, 92.2, 92.2, 92.2, 92.2, 92.2, 92.2, 92.2, 46.1]

# sampling fluctuations for 100, 200 and 300 um Si
# generation_sampling = [0.55, 2.2, 6.6]
# noise in mips for 200 um
generation_noise_sigma = [0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131, 0.131,
                          0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
