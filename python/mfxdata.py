"""
Utility class to read MinFluX data from Abberior microscopes
specpy class is from Imspector installation (C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\)
Author: Antonio Politi, MPINAT, 07.2022
"""
# TODO save as mat file for clustering (this should move to python also)

import zarr
import os
from specpy import *
import numpy as np
from scipy.spatial import transform
from scipy.interpolate import interp1d
import math
from matplotlib import pyplot as plt


class MfxData:
    MAX_TDIFF_REF = 10  # time difference between final record of ref (beads) and mfx recording in sec.
    # mfx recording is cropped accordingly
    ref_all = {}
    mfx_all = {}
    valid_ref_beads = {}  # Beads recording that fulfill minimal requirements
    maxtime_ref_beads = {}  # Consensus max time between reference beads recording
    TRANS = 'translate'
    ROT = 'rotate'

    def __init__(self, msrfile):
        if not os.path.exists(msrfile):
            raise OSError(msrfile + " file does not exist!")
        self.msrfile = msrfile
        self.filepath = os.path.dirname(msrfile)
        self.filename = os.path.splitext(os.path.basename(self.msrfile))[0]
        self.outdir = os.path.join(self.filepath, self.filename)

    def zarr_export(self, outdir=None, force=False):
        """
        Export mfx msr file to zarr files if data does not exist yet.
        :param outdir: path string. If None data is written to directory of msr file/msr_filename
        :param force: If zarr data exists force export
        :return: None
        """
        if outdir is not None:
            self.outdir = os.path.join(self.filepath, self.filename)
        afile = specpy.File(self.msrfile, File.Read)
        mf_data_sets = afile.minflux_data_sets()
        self.zarrdir = {mf['label']: os.path.join(self.outdir, mf['label']) for mf in mf_data_sets}
        # Sort assuming that order of acquisition is reflected in name
        self.zarrdir = dict(sorted(self.zarrdir.items()))
        for mf in mf_data_sets:
            if os.path.exists(self.zarrdir[mf['label']]) and not force:
                continue
            else:
                afile.unpack(mf['sid'], self.zarrdir[mf['label']])

    def zarr_import(self, force=False):
        """
        Import zarr files located in self.zarrdir if no data has been loaded yet
        :param force: force import
        """
        if len(self.ref_all) == 0 or force:
            for label, adir in self.zarrdir.items():
                zarr_data = zarr.open(store=os.path.join(adir, 'zarr'), mode='r')
                self.ref_all[label] = zarr_data.grd.mbm
                self.mfx_all[label] = zarr_data.mfx

    def get_valid_ref(self, force=False):
        """
        Make consistency check on ref_all data. And check that the beads are correct and remains attached throughout the experiment.
        :return: valid_ref_beads: A dictionary of label/experiment_id and bead.
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

        return valid_ref_beads

    def get_ref_vector(self):
        self.get_valid_ref()
        pos_array = []
        time_vector = []
        time_std_vector = []

        # Check for length of recording and consistency in difference between measurements
        for label, beads_id in self.valid_ref_beads.items():
            n_el = [len(self.ref_all[label][b]) for b in beads_id]
            if n_el.count(n_el[0]) != len(n_el):
                n_el = min(n_el)
            else:
                n_el = n_el[0]
            # Average time through the measurements of each round
            time_stack = np.stack([self.ref_all[label][b]['tim'][:n_el] for b in beads_id], axis=1)
            time_mean = np.mean(time_stack, axis=1)
            self.maxtime_ref_beads[label] = time_mean[-1]
            time_vector.append(time_mean)
            time_std_vector.append(np.std(time_stack, axis=1))
            pos_array.append(np.stack([self.ref_all[label][b]['pos'][:n_el] for b in beads_id], axis=1))

        # TODO This is not consistent in case there is a variation in order of label and idx
        for i in range(1, len(time_vector)):
            time_vector[i] = time_vector[i] + time_vector[i - 1][-1]

        time_vector = np.concatenate(time_vector)
        time_std_vector = np.concatenate(time_std_vector)
        pos_array = np.concatenate(pos_array)
        return [time_vector, pos_array, time_std_vector]

    def compute_ref_transform_error(self, pos_array):
        # Compute error in bead position now just std

        std = [np.std(pos_array[:, :, ax], axis=0) for ax in [0, 1, 2]]
        return {'ste': np.mean(np.linalg.norm(std, axis=1)) / math.sqrt(len(std)),
                'std': np.mean(np.linalg.norm(std, axis=1))}

    def compute_ref_transform(self):
        # Get vector of positions and time
        [time_vector, pos_array, time_std_vector] = self.get_ref_vector()

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

        out_dict = {}
        out_dict[self.TRANS] = translate
        out_dict[self.ROT] = rotate
        return out_dict

    def apply_ref_transform(self, translate, rotate, pos_array, time_vector):
        # Translate
        pos_array_reg = pos_array - translate(time_vector)
        # Rotate
        for idx in range(pos_array_reg.shape[0]):
            pos_array_reg[idx] = rotate(time_vector[idx]).apply(pos_array_reg[idx])
        return pos_array_reg

    def show_ref_transform(self, translate=None, rotate=None):

        if translate is None or rotate is None:
            transform = self.compute_ref_transform()
            translate = transform[self.TRANS]
            rotate = transform[self.ROT]

        [time_vector, pos_array, time_std_vector] = self.get_ref_vector()
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
        #

        transform = self.compute_ref_transform()
        lnc = {}
        tim = {}
        tid = {}
        loc = {}
        loc_reg = {}
        keys = list(self.mfx_all)
        for label, obj in self.mfx_all.items():
            # Trim on valid tracks
            lnc[label] = obj['itr']['lnc'][obj['vld']]
            loc[label] = obj['itr']['loc'][obj['vld']]
            tim[label] = obj['tim'][obj['vld']]
            tid[label] = obj['tid'][obj['vld']]

            # further trim for time ref_beads
            tim_trim = tim[label] < self.maxtime_ref_beads[label]
            lnc[label] = lnc[label][tim_trim]
            loc[label] = loc[label][tim_trim]
            tim[label] = tim[label][tim_trim]
            tid[label] = tid[label][tim_trim]
            # Keep only last iteration
            lnc[label] = lnc[label][:, -1]
            loc[label] = loc[label][:, -1]

        # Add time
        addtime = 0
        for idx in range(1, len(keys)):
            addtime += self.maxtime_ref_beads[keys[idx - 1]]
            tim[keys[idx]] = tim[keys[idx]] + addtime

        for label in self.mfx_all:
            loc_reg[label] = self.apply_ref_transform(transform[self.TRANS], transform[self.ROT], lnc[label],
                                                      tim[label])
        return [tim, lnc, loc, loc_reg]


if __name__ == "__main__":
    mfx = MfxData("C:/Users/apoliti/Desktop/mfluxtest/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.msr")
    mfx.zarr_export()
    outtrans = mfx.align_to_ref()s
    transform = mfx.compute_ref_transform()
    mfx.show_ref_transform(transform[mfx.TRANS], transform[mfx.ROT])
