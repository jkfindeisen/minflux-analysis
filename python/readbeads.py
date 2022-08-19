import zarr
import os
import numpy as np
import matplotlib.pyplot as plt

afile = specpy.File(msrfile, File.Read)

# get minflux data sets IDs save with their labels
mf_data_sets = afile.minflux_data_sets()
for mf in mf_data_sets:
    direxp = os.path.join(outdir, mf['label'])
    afile.unpack(mf['sid'], direxp)

from tkinter.filedialog import askdirectory
zarr_directory = askdirectory(initialdir='C:/Data')
import napari
napari.view_image(imga).
beads = zarr.open(store = zarr_directory, mode='r')

for R in beads['grd/mbm']:
    pos = beads['grd/mbm/' + R]['pos']
    np.savetxt("C:/Users/apoliti/my_data" + R + ".csv", pos, delimiter=",")

