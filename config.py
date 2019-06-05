#############################
# File reader
#
baseline_correct_npoints = 200 # Set to False to remove correction


#############################
# Time residual scripts
#
# Channel BandWidths to be applied in digital filtering
signal_BW      = 500e6  # Hz
trigger_BW     = 300e6  # Hz
NEMO_BW        = 50e6   # Hz

# Peak finding thresholds
signal_thresh  = -0.03  # V
trigger_thresh = -0.2   # V
NEMO_thresh    = -0.015 # V

# Integration window for charge measurement (applied about first peak)
signal_window  = [-10, 20]   # ns
trigger_window = [-10, 30]   # ns
NEMO_window    = [-10, 150]  # ns
