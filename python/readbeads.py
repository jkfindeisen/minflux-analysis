import zarr
import os
import numpy as np
from tkinter.filedialog import askdirectory
zarr_directory = askdirectory(initialdir='C:/Data')

beads = zarr.open(store = zarr_directory, mode='r')

for R in beads['grd/mbm']:
    pos = beads['grd/mbm/' + R]['pos']
    np.savetxt("C:/Users/apoliti/my_data" + R + ".csv", pos, delimiter=",")

