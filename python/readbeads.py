import zarr
import os
import numpy as np
import matplotlib.pyplot as plt
plt
from tkinter.filedialog import askdirectory
zarr_directory = askdirectory(initialdir='C:/Data')
import napari
napari.view_image(imga).
beads = zarr.open(store = zarr_directory, mode='r')

for R in beads['grd/mbm']:
    pos = beads['grd/mbm/' + R]['pos']
    np.savetxt("C:/Users/apoliti/my_data" + R + ".csv", pos, delimiter=",")

