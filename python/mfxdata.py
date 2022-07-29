"""
Utility class to read MinFluX data from Abberior microscopes
specpy class is from Imspector installation (C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\)
Author: Antonio Politi, MPINAT, 07.2022
"""


import zarr
import os
from specpy import *
import numpy as np
from scipy.interpolate import *
from scipy.linalg import norm
from scipy.spatial.transform import Rotation as R
from matplotlib import pyplot as plt
class MfxData:
    MAX_TDIFF_REF = 10 # time difference between final record of ref (beads) and mfx recording
    ref_all = {}
    mfx_all = {}

    def __init__(self, msrfile):
        if not os.path.exists(msrfile):
            raise OSError(msrfile + " file does not exist!")
        self.msrfile = msrfile
        self.filepath = os.path.dirname(msrfile)
        self.filename = os.path.splitext(os.path.basename(self.msrfile))[0]
        self.outdir = os.path.join(self.filepath, self.filename)

    def zarr_export(self, outdir = None, force = False):
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

    def zarr_import(self, force = False):
        if len(self.ref_all) == 0 or force:
            for label, adir in self.zarrdir.items():
                zarr_data = zarr.open(store = os.path.join(adir, 'zarr'), mode='r')
                self.ref_all[label] = zarr_data.grd.mbm
                self.mfx_all[label] = zarr_data.mfx

    def valid_ref(self):
        self.zarr_import()
        valid_ref_beads = {}
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
        return valid_ref_beads

    def compute_ref_transform(self):

        vld_ref = self.valid_ref()
        # Concatenate labels
        pos_array = []
        time_array = []
        #vld_ref = {'P1': ['R9', 'R10'], 'P2': ['R9', 'R10']}
        for label, beads_id in vld_ref.items():
            # check length of recording
            n_el = [len(self.ref_all[label][b]) for b in beads_id]
            if n_el.count(n_el[0]) != len(n_el):
                n_el = min(n_el)
            else:
                n_el = n_el[0]
            # create a column wise matrix for the time and position.
            time_array.append(np.mean(np.stack([self.ref_all[label][b]['tim'][:n_el] for b in beads_id], axis=1), axis=1))
            pos_array.append(np.stack([self.ref_all[label][b]['pos'][:n_el] for b in beads_id], axis=1))
        # Add time of last recording of beads
        time_array[1] = time_array[1] + time_array[0][-1]
        time_array = np.concatenate(time_array)
        pos_array = np.concatenate(pos_array)


        # TODO: if >=3 valid beads can also compute rotation matrix using Kabsch_algorithm scipy.spatial.transform.Rotation.align_vectors
        # Use centroids compute transform reduces std position for each bead

        centroids = np.mean(pos_array, axis=1)
        dpos = centroids - centroids[0]
        dpos_array = np.reshape(np.repeat(dpos, pos_array.shape[1], axis = 1 ), pos_array.shape, order = 'F')
        pos_array_translate = (pos_array - dpos_array)
        centroids_i = np.mean(pos_array_translate, axis=1)
        stdpos = [np.std(pos_array, axis = 0), np.std(pos_array_translate, axis=0)]

        # Try translation and rotation
        centroids = np.mean(pos_array, axis=1)
        dpos_array = np.reshape(np.repeat(centroids, pos_array.shape[1], axis = 1 ), pos_array.shape, order = 'F')

        pos_array_translate_rotate = (pos_array - dpos_array)

        for idx in range(pos_array_translate_rotate.shape[0]):
            rot = R.align_vectors(pos_array_translate_rotate[0], pos_array_translate_rotate[idx])
            pos_array_translate_rotate[idx] = rot[0].apply(pos_array_translate_rotate[idx])

        stdpos = [np.std(pos_array, axis = 0), np.std(pos_array_translate, axis=0), np.std(pos_array_translate_rotate, axis=0)]
        fig, axs = plt.subplots(3, 3, sharex=True, sharey = True)
        labels =['x_pos', 'y_pos', 'z_pos']
        for ax_idx, ax  in enumerate(axs):
            # Loop through the reference objects
            for idx in range(pos_array.shape[1]):
                ax[0].plot(time_array, pos_array[:, idx, ax_idx] - pos_array[0, idx, ax_idx])
                ax[1].plot(time_array, pos_array_translate[:, idx, ax_idx] - np.mean(pos_array_translate[0, idx, ax_idx]), ls="--", label=vld_ref['P1'][idx])
                ax[2].plot(time_array, pos_array_translate_rotate[:, idx, ax_idx] - np.mean(pos_array_translate_rotate[0, idx, ax_idx]), ls="--", label=vld_ref['P1'][idx])
            ax[2].legend()

            ax[0].set_ylabel(labels[ax_idx])

        axs[0][0].set_title('Unregistered')
        axs[0][1].set_title('Registered translate')
        axs[0][2].set_title('Registered translate + rotate')


        plt.show()

    def align_to_ref(self):
        #

        valid_ref_beads = self.valid_ref()

        zarr_data = zarr.open(store = os.path.join(self.zarrdir['P1'], 'zarr'), mode='r')
        ref = zarr_data.grd.mbm
        mfx = zarr_data.mfx





if __name__ == "__main__":
    mfx = MfxData("C:/Users/apoliti/Desktop/mfluxtest/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.msr")
    mfx.zarr_export()
    print(mfx.valid_ref())
    mfx.compute_ref_transform()

