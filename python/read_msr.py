from specpy import *
import zarr
import os

outdir = 'C:/Users/apoliti/Desktop/mfluxtest'
msrfile = "C:/Users/apoliti/Desktop/mfluxtest/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.msr"
basename = os.path.splitext(os.path.basename(msrfile))[0]
outdir = os.path.join(outdir, basename)


if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

afile = specpy.File(msrfile, File.Read)

# get minflux data sets IDs save with their labels
mf_data_sets = afile.minflux_data_sets()
for mf in mf_data_sets:
    direxp = os.path.join(outdir, mf['label'])
    afile.unpack(mf['sid'], direxp)

beads = zarr.open(store = zarr_directory, mode='r')

#zarr.f
# afile.read(14) stack corresponding to image.






