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

class MfxData:
    MAX_TDIFF_REF = 10    # time difference between final record of ref (beads) and mfx recording in sec.
    CORRECT_Z_POSITION_FACTOR = 0.7 # Correction due to oil water RI difference
    CORRECT_Z_POSITION_FACTOR_REF = 0.7

    ref_all = {}          # stores ref beads
    mfx_all = {}          # stores mfx measurements all washes
    valid_ref_beads = {}  # Beads recording that fulfill minimal requirements
    maxtime_ref_beads = {}  # Consensus max time between reference beads recording. Mfx recording is cropped accordingly
    msrfile_path = None
    msrfile_name = None
    outdir = None
    zarrdir = None
    TRANS = 'translate'
    ROT = 'rotate'

    def __init__(self, file_path):
        if not os.path.exists(file_path):
            raise OSError(file_path + " file does not exist!")
        self.msrfile_path = file_path

        self.msrfile_name = os.path.splitext(os.path.basename(self.msrfile_path))[0]
        # These are default values that can be overwritten

        self.outdir = self.get_outdir(self.msrfile_path)
        self.zarrdir = self.get_zarrdir(self.outdir, self.msrfile_path)


    @staticmethod
    def get_outdir(msrfile):
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
        zarrdir = dict(sorted(zarrdir.items()))
        return zarrdir

    def zarr_export(self):
        """
        Export mfx msr file to zarr files if data does not exist yet.
        :param zarrdir: A dictionary containing per label the path to export the msr file
        :param msrfile: path of msrfile
        :return: None
        """
        afile = specpy.File(self.msrfile_path, File.Read)
        mf_data_sets = afile.minflux_data_sets()
        for mf in mf_data_sets:
            afile.unpack(mf['sid'], self.zarrdir[mf['label']])

    def zarr_import(self, force=False):
        """
        Import zarr files located in self.zarrdir if no data has been loaded yet
        :param zarrdir: A dictionary containing per label the path to import the msr file
        :param msrfile: path of msrfile
        :param force: force import
        """

        # check that zarr files exist and eventually export
        for label, adir in self.zarrdir.items():
            if not os.path.exists(adir):
                self.zarr_export()

        if len(self.ref_all) == 0 or len(self.mfx_all) == 0 or force:
            for label, adir in self.zarrdir.items():
                zarr_data = zarr.open(store=os.path.join(adir, 'zarr'), mode='r')
                self.ref_all[label] = zarr_data.grd.mbm
                self.mfx_all[label] = zarr_data.mfx



    def set_valid_ref(self, force=False):
        """
        Make consistency check on ref_all, e.g. long enough recording.
        internal variable will be set in this function
        :param force: force import
        """
        # TODO When unequal number of references find consensus
        self.zarr_import()
        valid_ref_beads = {}
        if force or len(self.valid_ref_beads) > 0:
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
        """
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
            time_stack = np.stack([self.ref_all[label][b][col.TIM][:n_el] for b in self.valid_ref_beads], axis=1)
            time_mean = np.mean(time_stack, axis=1)
            self.maxtime_ref_beads[label] = time_mean[-1] # CONSISTENCY!
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
        # Compute error in bead position now just std

        std = [np.std(pos_array[:, :, ax], axis=0) for ax in [0, 1, 2]]
        return {'ste': np.mean(np.linalg.norm(std, axis=1)) / math.sqrt(len(std)),
                'std': np.mean(np.linalg.norm(std, axis=1))}

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
        titles = ['Unregistered, std (nm) %.2f, se (nm) %.2f\n', 'Translate, std (nm) %.2f, se (nm) %.2f\n', 'Translate + rotate, std (nm) %.2f, se (nm) %.2f\n']
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
            axs[0][iobj].set_title(titles[iobj] % (ref_error[iobj]['std'] * math.pow(10, 9),
                                                   ref_error[iobj]['ste'] * math.pow(10, 9)), fontsize=8)
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
            # further trim for time ref_beads
            tim_trim = tim < self.maxtime_ref_beads[label]
            lnc = lnc[tim_trim]
            loc = loc[tim_trim]
            tim = tim[tim_trim]
            tid = tid[tim_trim]
            vld = vld[tim_trim]
            # Keep only last iteration
            lnc = lnc[:, -1]
            loc = loc[:, -1]
            lnc[:, 2] = lnc[:, 2]*self.CORRECT_Z_POSITION_FACTOR
            loc[:, 2] = loc[:, 2]*self.CORRECT_Z_POSITION_FACTOR
            tim_tid_mean = self.get_meantime_trackid(tid, tim)
            out_dic[label] = {col.TIM: tim, col.TIM_TID_MEAN: tim_tid_mean, col.TID: tid, col.LNC: lnc,
                              col.LOC: loc, col.VLD: vld}
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

    def export_vtu(self, out_dict, lcoord, file_path):
        # export to vtk format for view in paraview, ideally also clustering etc.

        pos_concat = np.concatenate([out_dict[d][lcoord] for d in out_dict])
        tid_concat = np.concatenate([out_dict[d][col.TID] for d in out_dict])
        tim_concat = np.concatenate([out_dict[d][col.TIM] for d in out_dict])
        # concatenate positions
        x = np.ascontiguousarray(pos_concat[:, 0], dtype = np.float64)
        y = np.ascontiguousarray(pos_concat[:, 1], dtype = np.float64)
        z = np.ascontiguousarray(pos_concat[:, 2], dtype = np.float64)
        r = np.random.uniform(0, 1, x.shape[0]) # This needs to be computed
        keys = list(out_dict.keys())
        p = np.repeat(range(1, len(keys)+1), [out_dict[k][lcoord].shape[0] for k in keys])
        hl.pointsToVTK(file_path, x, y, z, data={'P': p, col.TID: tid_concat, col.TIM:  tim_concat,
                                                 'radius': r})

if __name__ == "__main__":
    mfx = MfxData("C:/Users/apoliti/Desktop/mfluxtest/data/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.msr")
    # different default output
    mfx.outdir = os.path.join("C:/Users/apoliti/Desktop/mfluxtest/analysis", mfx.msrfile_name)
    mfx.zarrdir = mfx.get_zarrdir(mfx.outdir, mfx.msrfile_path)

    mfx.zarr_import()
    mfx.set_valid_ref()

    regis = mfx.get_ref_transform()
    out_dic = mfx.align_to_ref()
    mfx.export_vtu(out_dic, col.LOC, "C:/Users/apoliti/Desktop/TMP/points_loc")
    mfx.export_vtu(out_dic, col.LTR, "C:/Users/apoliti/Desktop/TMP/points_ltr")
    mfx.export_vtu(out_dic, col.LRE, "C:/Users/apoliti/Desktop/TMP/points_lre")

    mfx.export_numpy(out_dic)
    #mfx.export_mat(out_dic)
    #mfx.export_ref_mat()
    #mfx.show_ref_transform(regis[mfx.TRANS], regis[mfx.ROT], save=True, show=False)
