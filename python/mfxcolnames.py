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
TID2 = 'tid2'           # Track ID given to the data internally. This can account also for splits events
STD = 'std'             # standard deviation through all dimensions sqrt(1/3(std_x^2 + std_y^2 + std_z^2))
SE = 'se'              # standard error through all dimensions   std/sqt(NL)
NL = 'nl'               # number of localizations
STD_XYZ = 'std_xyz'         # standard deviation in [std_x, std_y, std_z]
SE_XYZ = 'se_xyz'           # standard error in X