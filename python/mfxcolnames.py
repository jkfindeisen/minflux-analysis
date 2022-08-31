# Original columns (not complete)

LNC = 'lnc'             # Coordinates not registered
LOC = 'loc'             # From Imspector registered localization
VLD = 'vld'             # Valid recording
TID = 'tid'             # track ID
TIM = 'tim'             # time of recording
ITR = 'itr'             # iteration object
POS = 'pos'             # coordinate for reference beads

# Additional columns added to the data

LTR = 'ltr'             # Coordinates registered by translation using reference beads throughout the sample
LRE = 'lre'             # Coordinates registered by translation + rotation
TIM_TID_MEAN = 'tim_tid_mean' # Average acquisition time for all localizations within a track
CLS_TRACK = 'cls_track' # Cluster column for each tid
CLS_MEAS = 'cls_meas'   # Cluster column for each measurement
CLS_ALL = 'cls_all'     # Cluster column for joint measurement
