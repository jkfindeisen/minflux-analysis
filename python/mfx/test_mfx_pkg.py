from mfx.mfxdata import MfxData
from mfx.processlocalizations import ProcessLocalizations
import mfx.mfxcolnames as col
import time

# t0 = time.time()
# loc_dir = 'C:/Users/apoliti/Desktop/mflux_zarr_tmp_storage/analysis/Multiwash/VGLUT1_VGLUT1/'
# glob_dir = 'Z:/siva_minflux/analysis/Multiwash/VGLUT1_VGLUT1/'
# mfx = MfxData('Z:/siva_minflux/data/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01_First.msr',
#               outdir_main=glob_dir, zarr_dir_main=loc_dir)
#
# # Example of workaround for not merged msr file
# mfx.zarrdir = {'220811_VGLUT1_P1': 'C:/Users/apoliti/Desktop/mflux_zarr_tmp_storage/analysis/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01/220811_VGLUT1_P1',
#                '220811_VGLUT1_P2': 'C:/Users/apoliti/Desktop/mflux_zarr_tmp_storage/analysis/Multiwash/VGLUT1_VGLUT1/220811_VGLUT1_ROI01/220811_VGLUT1_P2'}
# mfx.zarr_import()
# mfx.set_valid_ref()
# regis = mfx.get_ref_transform()
# out_dic = mfx.align_to_ref()
# mfx.show_ref_transform(regis[mfx.TRANS], rotate=None, save=False, show=True)


pl = ProcessLocalizations(
    'Z:/siva_minflux/analysis//Multiwash/ZnT3_Syp/220309_ZnT3_Syp_ROI02/220309_ZnT3_Syp_ROI02.npy')


pl.log_parameters()
pl.trim_min_localizations()
pl.cluster_tid(method=pl.CLS_METHOD_BAYES_GMM)
pl.get_split_events()

#pl.DBCLUSTER_EPS_MEAS = 4e-9
#pl.cluster_meas(method=pl.CLS_METHOD_DBSCAN)
#pl.cluster_all(method=pl.CLS_METHOD_DBSCAN)
#pl.cluster_all_intersect(col.CLS_ALL)
#summary_dict = pl.summary_per_tid2()
#pl.export_vtu(summary_dict, col.LTR, file_path =pl.file_path_no_ext() + 'tmp', unit=col.UNIT_UM)
#summary_dict_flt = pl.filter_summary_loc(summary_dict,  filter_cluster={'name': col.CLS_MERGED_MEAS, 'size': 3})

#pl.export_vtu(summary_dict_flt, col.LTR, file_path =pl.file_path_no_ext() + '_tmp_flt_merged_meas_size3',
#              unit=col.UNIT_UM)

#pl.DBCLUSTER_EPS_MERGED_MEAS = 3e-8
#pl.DBCLUSTER_EPS_MERGED_ALL = 3e-8
#pl.log_parameters()
#summary_dict = pl.summary_per_tid2()
#pl.export_csv(summary_dict, file_path=pl.file_path_no_ext())
#pl.export_vtu(in_dict=summary_dict, coord=col.LTR, file_path=pl.file_path_no_ext())

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
