"""
Utility class to read MinFluX data from Abberior microscopes
specpy class is from Imspector installation (C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\)
Author: Antonio Politi, MPINAT, 07.2022
"""
# TODO save as mat file for clustering (this should move to python also)
import scipy.io
import zarr
import os
from specpy import *
import numpy as np
from scipy.spatial import transform
from scipy.interpolate import interp1d
import math
from matplotlib import pyplot as plt


class MfxData:
    MAX_TDIFF_REF = 10    # time difference between final record of ref (beads) and mfx recording in sec.
    ref_all = {}          # stores ref beads
    mfx_all = {}          # stores mfx measurements all washes
    valid_ref_beads = {}  # Beads recording that fulfill minimal requirements
    maxtime_ref_beads = {}  # Consensus max time between reference beads recording. Mfx recording is cropped accordingly


    TRANS = 'translate'
    ROT = 'rotate'

    def __init__(self, msrfile):
        if not os.path.exists(msrfile):
            raise OSError(msrfile + " file does not exist!")
        self.msrfile = msrfile

        self.filepath = os.path.dirname(msrfile)
        self.filename = os.path.splitext(os.path.basename(self.msrfile))[0]
        # These are default values that can be overwritten

        self.outdir = self.get_outdir(self.msrfile)
        self.zarrdir = self.get_zarrdir(self.outdir, self.msrfile)


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
        afile = specpy.File(self.msrfile, File.Read)
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
        self.zarr_import()
        valid_ref_beads = {}
        if len(self.valid_ref_beads) == 0 or force:
            for label in self.ref_all:
                ref = self.ref_all[label]
                mfx = self.mfx_all[label]
                valid_tmax = max(mfx['tim'][mfx['vld']])
                vld = []
                for key, r in ref.items():
                    if len(r['tim']) == 0:
                        continue
                    if max(r['tim']) > valid_tmax:
                        vld.append(key)
                        continue
                    if abs(max(r['tim']) - valid_tmax) <= self.MAX_TDIFF_REF:
                        vld.append(key)
                valid_ref_beads[label] = vld
            self.valid_ref_beads = valid_ref_beads

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
        for label, beads_id in self.valid_ref_beads.items():
            # Have same length
            n_el = [len(self.ref_all[label][b]) for b in beads_id]
            if n_el.count(n_el[0]) != len(n_el):
                n_el = min(n_el)
            else:
                n_el = n_el[0]
            # Average time through the measurements of each round
            time_stack = np.stack([self.ref_all[label][b]['tim'][:n_el] for b in beads_id], axis=1)
            time_mean = np.mean(time_stack, axis=1)
            self.maxtime_ref_beads[label] = time_mean[-1] # CONSISTENCY!
            time_vector.append(time_mean)
            time_std_vector.append(np.std(time_stack, axis=1))
            pos_array.append(np.stack([self.ref_all[label][b]['pos'][:n_el] for b in beads_id], axis=1))

        # TODO This is not consistent in case there is a variation in order of label and idx ??
        for i in range(1, len(time_vector)):
            time_vector[i] = time_vector[i] + time_vector[i - 1][-1]

        time_vector = np.concatenate(time_vector)
        time_std_vector = np.concatenate(time_std_vector)
        pos_array = np.concatenate(pos_array)
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
        pos_array_reg = self.apply_ref_rotate(rotate, pos_array_reg, time_vector)
        return pos_array_reg

    def show_ref_transform(self, translate, rotate):
        [time_vector, pos_array, time_std_vector] = self.get_ref()
        pos_array_translate_rotate = np.zeros_like(pos_array) # this is important otherwise pas by reference
        for idx in range(0, pos_array.shape[1]):
            pos_array_translate_rotate[:, idx] = self.apply_ref_transform(translate, rotate, pos_array[:, idx],
                                                                          time_vector)

        fig, axs = plt.subplots(3, 2, sharex=True, sharey=True)
        labels = ['x_pos', 'y_pos', 'z_pos']
        for ax_idx, ax in enumerate(axs):
            # Loop through the reference objects
            for idx in range(pos_array.shape[1]):
                ax[0].plot(time_vector, pos_array[:, idx, ax_idx] - pos_array[0, idx, ax_idx])
                ax[1].plot(time_vector, pos_array_translate_rotate[:, idx, ax_idx] -
                           np.mean(pos_array_translate_rotate[0, idx, ax_idx]),
                           ls="--", label=self.valid_ref_beads['P1'][idx])

            ax[1].legend()
            ax[0].set_ylabel(labels[ax_idx])
        ref_error = self.compute_ref_transform_error(pos_array)
        ref_error_translate_rotate = self.compute_ref_transform_error(pos_array_translate_rotate)

        axs[0][0].set_title('Unregistered, std (nm) %.2f, se (nm) %.2f' % (ref_error['std'] * math.pow(10, 9),
                                                                           ref_error['ste'] * math.pow(10, 9)))
        axs[0][1].set_title('Reg. translate + rotate, std (nm) %.2f, se (nm) %.2f' % (
            ref_error_translate_rotate['std'] * math.pow(10, 9),
            ref_error_translate_rotate['ste'] * math.pow(10, 9)))
        plt.show()

    def align_to_ref(self):
        register = self.get_ref_transform()
        out_dic = {}
        keys = list(self.mfx_all)
        for label, obj in self.mfx_all.items():
            # Trim on valid tracks
            lnc = obj['itr']['lnc'][obj['vld']]
            loc = obj['itr']['loc'][obj['vld']]
            tim = obj['tim'][obj['vld']]
            tid = obj['tid'][obj['vld']]

            # further trim for time ref_beads
            tim_trim = tim < self.maxtime_ref_beads[label]
            lnc = lnc[tim_trim]
            loc = loc[tim_trim]
            tim = tim[tim_trim]
            tid = tid[tim_trim]
            # Keep only last iteration
            lnc = lnc[:, -1]
            loc = loc[:, -1]
            out_dic[label] = {'tim': tim, 'tid': tid, 'lnc': lnc, 'loc': loc}


        # Add time
        add_time = 0
        for idx in range(1, len(keys)):
            add_time += self.maxtime_ref_beads[keys[idx - 1]]
            out_dic[keys[idx]]['tim'] = out_dic[keys[idx]]['tim'] + add_time

        # register
        for label in self.mfx_all:
            out_dic[label]['lre'] = self.apply_ref_transform(register[self.TRANS], register[self.ROT],
                                                             out_dic[label]['lnc'], out_dic[label]['tim'])
        #for label in self.mfx_all:
        return out_dic

    def export_mat(self, out_dict):
        scipy.io.savemat(os.path.join(self.outdir, self.filename + ".mat"), out_dict)



if __name__ == "__main__":
    mfx = MfxData("C:/Users/apoliti/Desktop/mfluxtest/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.msr")
    mfx.zarr_import()

    regis = mfx.get_ref_transform()
    out_dic = mfx.align_to_ref()
    mfx.export_mat(out_dic)
    mfx.show_ref_transform(regis[mfx.TRANS], regis[mfx.ROT])
