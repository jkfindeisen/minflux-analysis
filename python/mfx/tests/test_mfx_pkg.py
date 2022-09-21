from ..mfx.mfxdata import MfxData
import time

t0 = time.time()
loc_dir = 'C:/Users/apoliti/Desktop/mflux_zarr_tmp_storage/analysis/Multiwash/VGLUT1_VGLUT1/'
glob_dir = 'Z:/siva_minflux/analysis/Multiwash/VGLUT1_VGLUT1/'
mfx = MfxData('Z:/siva_minflux/data/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01_First.msr',
              outdir_main=glob_dir, zarr_dir_main=loc_dir)

# Example of workaround for not merged msr file
mfx.zarrdir = {'220811_VGLUT1_P1': 'C:/Users/apoliti/Desktop/mflux_zarr_tmp_storage/analysis/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01/220811_VGLUT1_P1',
               '220811_VGLUT1_P2': 'C:/Users/apoliti/Desktop/mflux_zarr_tmp_storage/analysis/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01/220811_VGLUT1_P2'}
mfx.zarr_import()
mfx.set_valid_ref()
regis = mfx.get_ref_transform()
out_dic = mfx.align_to_ref()
mfx.show_ref_transform(regis[mfx.TRANS], rotate=None, save=False, show=True)

#mfx1.MAX_TDIFF_REF = 10
#mfx2 = MfxData('Z:/siva_minflux/data/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01_Second.msr',
#                outdir_main=glob_dir, zarr_dir_main=loc_dir)
#mfx2.MAX_TDIFF_REF = 10

#mfx1.zarr_import()
#mfx2.zarr_import()
#mfx1.set_valid_ref()
#mfx2.set_valid_ref()

# Match beads stump code

#[time_vector, pos_array1, time_std_vector] = mfx1.get_ref()
#[time_vector, pos_array2, time_std_vector] = mfx2.get_ref()
#mat_dd = distance_matrix(pos_array1[0], pos_array2[0])*1e9
#match_beads = list()
#max_beads_dist = 150
#for idx in range(0, len(mfx1.valid_ref_beads)):
#    if any(mat_dd[idx] < max_beads_dist):
#        match_beads.append([mfx1.valid_ref_beads[idx],
#                            mfx2.valid_ref_beads[np.where(mat_dd[idx] < max_beads_dist)[0][0]]])

#

#print("Time elapsed: ",  time.time() - t0)
#mfx.export_vtu(out_dic, col.LOC, "C:/Users/apoliti/Desktop/TMP/points_loc")
#mfx.export_vtu(out_dic, col.LTR, "C:/Users/apoliti/Desktop/TMP/points_ltr")
#mfx.export_vtu(out_dic, col.LRE, "C:/Users/apoliti/Desktop/TMP/points_lre")
#mfx.show_ref_transform(regis[mfx.TRANS], rotate=None, save=False, show=True)
#mfx.export_numpy(out_dic)
#mfx.export_mat(out_dic)
#mfx.export_ref_mat()
#
