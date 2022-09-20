"""
Utility class to read MinFluX data from Abberior microscopes and perform some post processing
The class can perform a registration using beads and combine different acquisition that have been subsequently performed (PAINT)
The registration transform is translation (recommended) or translation + rotation (!! Rotation gives wrong results if beads are colinear !!)
specpy class is from Imspector installation (C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\)
Author: Antonio Politi, MPINAT, 07.2022
"""
import scipy.io
import zarr
import os
from specpy import *
import numpy as np
from scipy.spatial import transform
from scipy.interpolate import interp1d
import math
from matplotlib import pyplot as plt
import warnings
from evtk import hl
import mfxcolnames as col
import re
import time
from scipy.spatial import distance_matrix


class MfxData:
    # class variable shared by all instances
    MAX_TDIFF_REF = 10    # time difference between final record of ref (beads) and mfx recording in sec.
    CORRECT_Z_POSITION_FACTOR = 0.7 # Correction due to oil water RI difference
    CORRECT_Z_POSITION_FACTOR_REF = 0.7
    TRANS = 'translate'
    ROT = 'rotate'
    LOAD_ZARR_IN_MEMORY = True

    def __init__(self, file_path, outdir_main=None, zarr_dir_main=None):
        # instance variables
        self.ref_all = {}          # stores ref beads
        self.mfx_all = {}          # stores mfx measurements all washes
        self.valid_ref_beads = {}  # Beads recording that fulfill minimal requirements
        self.maxtime_ref_beads = {}  # Consensus max time between reference beads recording. Mfx recording is cropped accordingly
        self.mintime_ref_beads = {}
        if not os.path.exists(file_path):
            raise OSError(file_path + " file does not exist!")
        self.msrfile_path = file_path

        self.msrfile_name = os.path.splitext(os.path.basename(self.msrfile_path))[0]
        # These are default values that can be overwritten

        self.outdir = self.get_outdir(self.msrfile_path, outdir_main)
        # Local directory performs 10x faster than network connection/ Staging on a local directory may be better
        if zarr_dir_main is None:
            zarr_dir_main = self.outdir
        else:
            zarr_dir_main = self.get_outdir(self.msrfile_path, zarr_dir_main)
        self.zarrdir = self.get_zarrdir(outdir=zarr_dir_main, msrfile=self.msrfile_path)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)


    @staticmethod
    def get_outdir(msrfile, filepath = None):
        if filepath == None:
            filepath = os.path.dirname(msrfile)
        filename = os.path.splitext(os.path.basename(msrfile))[0]
        # These are default values that can be overwritten
        return os.path.join(filepath, filename)

    @staticmethod
    def get_zarrdir(outdir, msrfile):
        afile = specpy.File(msrfile, File.Read)
        mf_data_sets = afile.minflux_data_sets()
        zarrdir = {mf['label']: os.path.join(outdir, mf['label']) for mf in mf_data_sets}
        # Sort assuming that order of acquisition is reflected in name
        zarrdir = dict(sorted(zarrdir.items(),  key=lambda x: int(re.search(r'\d+$', x[0]).group())))
        return zarrdir

    def __zarr_export(self):
        """
        Exports mfx msr file to zarr files if data does not exist yet.
        :return: None
        """
        afile = specpy.File(self.msrfile_path, File.Read)
        mf_data_sets = afile.minflux_data_sets()
        for mf in mf_data_sets:
            afile.unpack(mf['sid'], self.zarrdir[mf['label']])

    def zarr_import(self, force=False):
        """
        Import zarr files located in self.zarrdir if no data has been loaded yet. Eventually create zarr data structure
        :param force: force import
        """

        # check that zarr files exist and eventually export
        for label, adir in self.zarrdir.items():
            if not os.path.exists(adir):
                print('Create Zarr data structure')
                self.__zarr_export()

        if len(self.ref_all) == 0 or len(self.mfx_all) == 0 or force:
            for label, adir in self.zarrdir.items():
                dir_store1 = zarr.DirectoryStore(os.path.join(adir, 'zarr'))
                # dirr_store1.rename(src_dir=..., dst_dir=...) to change name of directory
                cache1 = zarr.LRUStoreCache(store=dir_store1, max_size=2**28)
                if self.LOAD_ZARR_IN_MEMORY:
                    # Does not seem to make a difference in speed
                    self.ref_all[label] = zarr.load(store=cache1, path='grd/mbm/')
                    self.mfx_all[label] = zarr.load(store=cache1, path='mfx')
                else:
                    self.ref_all[label] = zarr.open(store=cache1, mode='r',  path='grd/mbm/')
                    self.mfx_all[label] = zarr.open(store=cache1, mode='r', path='mfx')

    def set_valid_ref(self, force=False):
        """
        Make consistency check on ref_all, e.g. long enough recording.
        internal variable will be set in this function
        :param force: force import
        """
        # TODO: Create a match beads functionality if names are different
        valid_ref_beads = {}
        if not force and len(self.valid_ref_beads) > 0:
            return
        for label in self.ref_all:
            ref = self.ref_all[label]
            mfx = self.mfx_all[label]
            valid_tmax = max(mfx[col.TIM][mfx[col.VLD]])
            vld = []
            for key, r in ref.items():
                if len(r[col.TIM]) == 0:
                    continue
                if max(r[col.TIM]) > valid_tmax:
                    vld.append(key)
                    continue
                if abs(max(r[col.TIM]) - valid_tmax) <= self.MAX_TDIFF_REF:
                    vld.append(key)
            valid_ref_beads[label] = vld
        # Match the beads in the different washes
        keys = list(valid_ref_beads)

        vld_ref_all = valid_ref_beads[keys[0]]
        for i in range(1, len(keys)):
            vld_ref_all = [x for x in vld_ref_all if x in valid_ref_beads[keys[i]]]
        self.valid_ref_beads = vld_ref_all

    def get_ref(self):
        """
        Extract the positions of valid beads, concatenate those for all directories
        :return: a list with time_vector, pos_array of all beads, std in time in case there is a jump in the recording
        pos_array[time, bead, axes]
        """
        # TODO: Add dynamic time warping to match correctly the different recordings from all beads.
        self.set_valid_ref()
        pos_array = []
        time_vector = []
        time_std_vector = []

        # Check for length of recording and consistency in difference between measurements
        for label in self.ref_all:
            # Have same length
            n_el = [len(self.ref_all[label][b]) for b in self.valid_ref_beads]
            if n_el.count(n_el[0]) != len(n_el):
                n_el = min(n_el)
            else:
                n_el = n_el[0]
            # Average time through the measurements of each round
            # Rough find time shift between
            time_stack = np.stack([self.ref_all[label][b][col.TIM][:n_el] for b in self.valid_ref_beads], axis=1)
            time_mean = np.mean(time_stack, axis=1)
            self.maxtime_ref_beads[label] = time_mean[-1] # CONSISTENCY!
            self.mintime_ref_beads[label] = time_mean[0]
            time_vector.append(time_mean)
            time_std_vector.append(np.std(time_stack, axis=1))
            pos_array.append(np.stack([self.ref_all[label][b][col.POS][:n_el] for b in self.valid_ref_beads], axis=1))

        # TODO This is not consistent in case there is a variation in order of label and idx ??
        for i in range(1, len(time_vector)):
            time_vector[i] = time_vector[i] + time_vector[i - 1][-1]

        time_vector = np.concatenate(time_vector)
        time_std_vector = np.concatenate(time_std_vector)
        pos_array = np.concatenate(pos_array)
        pos_array[:, :, 2] = pos_array[:, :, 2]*self.CORRECT_Z_POSITION_FACTOR_REF
        return [time_vector, pos_array, time_std_vector]

    def compute_ref_transform_error(self, pos_array):
        """
        Some metrics to asses how stable the beads are
        :param pos_array:
        :return: the metric
        """
        # Compute error in bead position, this is the mean in all directions. This is the best fit of a symmetric gaussian
        # Alternatively one could compute the norm, which is the highest possible error!
        # The standard error is truly the error of the centroid ??

        std = [np.std(pos_array[:, :, ax], axis=0) for ax in range(0, 3)]
        std_m = np.mean(std, axis=1)
        # std[axis_dimension, bead_id]
        return {'se_xy': np.sqrt((std_m[0]**2 + std_m[1]**2)/2)/ math.sqrt(np.shape(std)[1]),
                'se_z': std_m[2]/math.sqrt(np.shape(std)[1]),
                'std_xy': np.sqrt((std_m[0]**2 + std_m[1]**2)/2),
                'std_z': std_m[2]}

    def get_ref_transform(self):
        """
        Compute translation and rotation
        :return:
        """
        # Vector of positions and time
        [time_vector, pos_array, time_std_vector] = self.get_ref()

        # Translation using the centroids
        centroids = np.mean(pos_array, axis=1)
        dpos_array = np.reshape(np.repeat(centroids, pos_array.shape[1], axis=1), pos_array.shape, order='F')
        pos_array_translate = (pos_array - dpos_array)
        # interpolate the translation
        translate = interp1d(time_vector, centroids, axis=0)

        if len(self.valid_ref_beads) < 3:
            warnings.warn("Less than 3 reference beads, rotation registration is not computed", category=UserWarning)
            rotate = None
        else:
            # Rotation
            rot = []
            for idx in range(pos_array_translate.shape[0]):
                [rot_loc, rmsd] = transform.Rotation.align_vectors(pos_array_translate[0], pos_array_translate[idx])
                rot.append(rot_loc)
            # concatenate rotations (newer version of scipy does it in one go)
            rot_vec = transform.Rotation.from_rotvec([r.as_rotvec() for r in rot])
            rotate = transform.Slerp(time_vector, rot_vec)

        out_transform = {self.TRANS: translate, self.ROT: rotate}
        return out_transform

    def apply_ref_translate(self, translate, pos_array, time_vector):
        return pos_array - translate(time_vector)

    def apply_ref_rotate(self, rotate, pos_array, time_vector):
        for idx in range(pos_array.shape[0]):
            pos_array[idx] = rotate(time_vector[idx]).apply(pos_array[idx])
        return pos_array

    def apply_ref_transform(self, translate, rotate, pos_array, time_vector):
        # apply translation followed by rotation
        pos_array_reg = self.apply_ref_translate(translate, pos_array, time_vector)
        pos_array_reg = self.apply_ref_rotate(rotate, pos_array_reg.copy(), time_vector)
        return pos_array_reg

    def show_ref_transform(self, translate, rotate, show=False, save=False):
        if not save and not show:
            return
        [time_vector, pos_array, time_std_vector] = self.get_ref()
        pos_array_translate = np.zeros_like(pos_array)
        pos_array_translate_rotate = np.zeros_like(pos_array) # this is important otherwise pas by reference

        for idx in range(0, pos_array.shape[1]):
            pos_array_translate[:, idx] = self.apply_ref_translate(translate, pos_array[:, idx],
                                                                   time_vector)
            if rotate is not None:
                pos_array_translate_rotate[:, idx] = self.apply_ref_transform(translate, rotate, pos_array[:, idx],
                                                                              time_vector)
        titles = ['Unregistered (nm), std_xy %.2f, std_z  %.2f\nse_xy %.2f, se_z %.2f\n',
                  'Translate, std_xy %.2f, std_z  %.2f\nse_xy %.2f, se_z %.2f\n',
                  'Translate + rotate, std_xy %.2f, std_z  %.2f\nse_xy %.2f, se_z %.2f\n']

        if rotate is None:
            obj_to_plot = [pos_array, pos_array_translate]
        else:
            obj_to_plot = [pos_array, pos_array_translate, pos_array_translate_rotate]

        fig, axs = plt.subplots(3, len(obj_to_plot), sharex=True, sharey=True, figsize=(10, 5))

        labels = ['x_pos', 'y_pos', 'z_pos']
        for ax_idx, ax in enumerate(axs):
            # Loop through the reference objects
            for idx in range(pos_array.shape[1]):
                for iobj in range(len(obj_to_plot)):
                    ax[iobj].plot(time_vector, obj_to_plot[iobj][:, idx, ax_idx] -
                                  obj_to_plot[iobj][0, idx, ax_idx], label=self.valid_ref_beads[idx])
            ax[0].legend(fontsize=8)
            ax[0].set_ylabel(labels[ax_idx])


        ref_error = [self.compute_ref_transform_error(obj) for obj in obj_to_plot]
        for iobj in range(len(obj_to_plot)):
            axs[0][iobj].set_title(titles[iobj] % (ref_error[iobj]['std_xy'] * math.pow(10, 9),
                                                   ref_error[iobj]['std_z'] * math.pow(10, 9),
                                                   ref_error[iobj]['se_xy'] * math.pow(10, 9),
                                                   ref_error[iobj]['se_z'] * math.pow(10, 9)
                                                   ), fontsize=8)
            axs[2][iobj].set_xlabel('time (sec)')
        if save:
            plt.savefig(os.path.join(self.outdir, self.msrfile_name + "_ref_drift.png"))
        if show:
            plt.show()

    def align_to_ref(self):
        register = self.get_ref_transform()
        out_dic = {}
        keys = list(self.mfx_all)
        for label, obj in self.mfx_all.items():
            # Trim on valid tracks
            print("Valid tracks %s %d/%d" % (label, len(np.unique(obj[col.TID][obj[col.VLD]])),
                                             len(np.unique(obj[col.TID]))))
            lnc_nv = obj[col.ITR][col.LNC][~obj[col.VLD]]
            lnc_nv = lnc_nv[:, -1]
            print("Localizations in invalid tracks %s %d/%d" % (label, len(lnc_nv[~np.isnan(lnc_nv[:, 0]), 0]),
                                                                len(lnc_nv[np.isnan(lnc_nv[:,0]), 0])))
            lnc = obj[col.ITR][col.LNC][obj[col.VLD]]
            loc = obj[col.ITR][col.LOC][obj[col.VLD]]
            tim = obj[col.TIM][obj[col.VLD]]
            tid = obj[col.TID][obj[col.VLD]]
            vld = obj[col.VLD][obj[col.VLD]]
            eco = obj[col.ITR][col.ECO][obj[col.VLD]]
            efo = obj[col.ITR][col.EFO][obj[col.VLD]]
            # further trim for time ref_beads
            tim_trim = (tim > self.mintime_ref_beads[label]) * (tim < self.maxtime_ref_beads[label])
            lnc = lnc[tim_trim]
            loc = loc[tim_trim]
            tim = tim[tim_trim]
            tid = tid[tim_trim]
            vld = vld[tim_trim]
            eco = eco[tim_trim]
            efo = efo[tim_trim]
            # Keep only last iteration
            lnc = lnc[:, -1]
            loc = loc[:, -1]
            eco = eco[:, -1]
            efo = efo[:, -1]
            lnc[:, 2] = lnc[:, 2]*self.CORRECT_Z_POSITION_FACTOR
            loc[:, 2] = loc[:, 2]*self.CORRECT_Z_POSITION_FACTOR
            tim_tid_mean = self.get_meantime_trackid(tid, tim)
            out_dic[label] = {col.TIM: tim, col.TIM_TID_MEAN: tim_tid_mean, col.TID: tid, col.LNC: lnc,
                              col.LOC: loc, col.VLD: vld, col.ECO: eco, col.EFO: efo}
        # Add time
        add_time = 0
        add_tid = 0
        for idx in range(1, len(keys)):
            add_tid += out_dic[keys[idx-1]][col.TID][-1]
            add_time += self.maxtime_ref_beads[keys[idx - 1]]
            out_dic[keys[idx]][col.TIM] = out_dic[keys[idx]][col.TIM] + add_time
            out_dic[keys[idx]][col.TIM_TID_MEAN] = out_dic[keys[idx]][col.TIM_TID_MEAN] + add_time
            out_dic[keys[idx]][col.TID] = out_dic[keys[idx]][col.TID] + add_tid

        # register
        for label in self.mfx_all:
            # Translate to center for convenience, translate each track with the same function
            out_dic[label][col.LNC] = self.apply_ref_translate(register[self.TRANS],
                                                               out_dic[label][col.LNC], out_dic[label][col.TIM][0])
            out_dic[label][col.LOC] = self.apply_ref_translate(register[self.TRANS],
                                                               out_dic[label][col.LOC], out_dic[label][col.TIM][0])
            # Translate over time
            out_dic[label][col.LTR] = self.apply_ref_translate(register[self.TRANS],
                                                               out_dic[label][col.LNC],
                                                               out_dic[label][col.TIM_TID_MEAN])
            if register[self.ROT] is None:
                out_dic[label][col.LRE] = np.zeros_like(out_dic[label][col.LNC])
            else:
                out_dic[label][col.LRE] = self.apply_ref_transform(register[self.TRANS], register[self.ROT],
                                                                   out_dic[label][col.LNC],
                                                                   out_dic[label][col.TIM_TID_MEAN])
        return out_dic

    def get_meantime_trackid(self, tid, tim):
        tout = np.copy(tim)
        for id in np.unique(tid):
            idx = np.where(tid == id)[0]
            tout[idx] = np.mean(tim[idx])
        return tout

    def export_numpy(self, out_dict):
        np.save(os.path.join(self.outdir, self.msrfile_name + ".npy"), out_dict)

    def export_mat(self, out_dict):
        scipy.io.savemat(os.path.join(self.outdir, self.msrfile_name + ".mat"), out_dict)

    def export_ref_mat(self):
        register = self.get_ref_transform()
        [time_vector, pos_array, time_std_vector] = self.get_ref()
        flat_pos = np.reshape(pos_array, (pos_array.shape[0]*pos_array.shape[1], pos_array.shape[2]))
        time_vector = np.repeat(time_vector, pos_array.shape[1])

        flat_pos_trans = self.apply_ref_translate(register[self.TRANS], flat_pos, time_vector)
        if register[self.ROT] is not None:
            flat_pos_reg = self.apply_ref_transform(register[self.TRANS], register[self.ROT],
                                                    flat_pos, time_vector)
        else:
            flat_pos_reg = np.zeros_like(flat_pos)
        # Center to 0 for rendering
        pos_array2 = pos_array.copy() - np.mean(pos_array, axis=1)[0]
        flat_pos2 = np.reshape(pos_array2, (pos_array2.shape[0]*pos_array2.shape[1], pos_array2.shape[2]))

        out_dict = {col.TIM: time_vector, col.LNC: flat_pos2, col.LTR: flat_pos_trans, col.LRE: flat_pos_reg}
        scipy.io.savemat(os.path.join(self.outdir, self.msrfile_name + "_ref.mat"), out_dict)



if __name__ == "__main__":
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
