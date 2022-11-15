import os
import glob
import warnings
from .mfxdata import MfxData
from .processlocalizations import ProcessLocalizations
from . import mfxcolnames as col
def gen_directories(main_dir, sub_dir, keys):
    """
    Combine different path to directories
    :param main_dir: a main directory
    :param sub_dir: a sub directory
    :param keys: a list of key values
    :return: out_dir_dict A dictionary where for each key it contains a path
    """
    dir_dict = {}
    for key in keys:
        dir_dict[key] = os.path.join(main_dir, sub_dir, key)
    return dir_dict


def gen_msrfiles(dir_dict):
    """
    Per each path/key value in dir_dict create list of msr files
    :param dir_dict: A dictionary of paths
    :return: npy_files: A dictionary of file arrays
    """
    msr_files = {}
    for key in dir_dict:
        files = glob.glob(os.path.join(dir_dict[key], '*.msr'))
        if len(files) == 0:
            warnings.warn("No msr files found in directory %s" % dir_dict[key])
        msr_files[key] = files
    return msr_files


def gen_npyfiles(dir_dict):
    """
    Per each path/key value in dir_dict create list of npy files
    :param dir_dict: A dictionary of paths
    :return: npy_files: A dictionary of file arrays
    """
    npy_files = {}
    for key in dir_dict:
        npy_files[key] = []

    for key in dir_dict:
        for (root, dirs, files) in os.walk(dir_dict[key]):
            if len(dirs) > 0:
                for adir in dirs:
                    npyfile = os.path.join(root, adir, adir + '.npy')
                    if os.path.exists(npyfile):
                        npy_files[key].append(npyfile)
    return npy_files


def process_mfx_npyfile(file_name, eps_range):
    """
    Example of workflow to trim and cluster localizations
    :param file_name:
    :param eps_range:
    """
    pl = ProcessLocalizations(file_name)
    pl.MIN_LOCALIZATIONS = 3    # 3 is the default value
    # Remove all tracks with less than MIN_LOCALIZATIONS
    pl.trim_min_localizations()
    # Split tracks if needed, remove outliers in each track
    pl.cluster_tid(method=pl.CLS_METHOD_BAYES_GMM)

    # find objects that belong to the same cluster according to distance 30 nm
    for eps in eps_range:
        pl.DBCLUSTER_EPS_MEAS = eps
        pl.DBCLUSTER_EPS_ALL = eps
        pl.DBCLUSTER_EPS_MERGED_ALL = eps
        pl.DBCLUSTER_EPS_MERGED_MEAS = eps
        pl.log_parameters()
        pl.cluster_meas(method=pl.CLS_METHOD_DBSCAN)
        pl.cluster_all(method=pl.CLS_METHOD_DBSCAN)
        pl.cluster_all_intersect(col.CLS_ALL)
        summary_dict = pl.summary_per_tid2()
        eps_txt = '%.0f' % (eps*1e9)
        pl.export_csv(summary_dict, file_path=pl.file_path_no_ext() + '_eps' + eps_txt)
        pl.export_vtu(in_dict=summary_dict, coord=col.LTR,
                      file_path=pl.file_path_no_ext() + '_eps' + eps_txt)


# Realign data funcion
def align_to_ref_beads(msr_files, out_dir, zarr_dir, invalid_msr, exclude_beads, time_diff_ref, wash_to_use=None):
    for wash_key in msr_files:
        if wash_to_use is not None and wash_key not in wash_to_use:
            continue
        for file_path in msr_files[wash_key]:

            file_name = os.path.basename(file_path)
            if file_name in invalid_msr:
                continue
            mfx = MfxData(file_path, outdir_main=out_dir[wash_key],
                                  zarr_dir_main=zarr_dir[wash_key])
            if file_name in time_diff_ref:
                mfx.MAX_TDIFF_REF = time_diff_ref[file_name]
            else:
                mfx.MAX_TDIFF_REF = 10
            mfx.zarr_import()
            mfx.set_valid_ref()
            print('\n' + mfx.msrfile_path)
            print(mfx.valid_ref_beads)
            if file_name in exclude_beads:
                mfx.valid_ref_beads = [x for x in mfx.valid_ref_beads if x not in exclude_beads[file_name]]
            print(mfx.valid_ref_beads)
            regis = mfx.get_ref_transform()
            mfx.show_ref_transform(translate=regis[mfx.TRANS], rotate=None, save=True, show=True)
            out_data_dict = mfx.align_to_ref()
            mfx.export_numpy(out_data_dict)
            mfx.export_ref_mat()


def preview_ref_beads(msr_files: dict, out_dir: dict, zarr_dir: dict, id_file: int, key: str,
                      time_diff_ref: float = 10, exclude_beads: list = None):
    """
    Preview ref beads and check their number and quality
    :param msr_files: Dictionary of msr files
    :param out_dir: Output directory path
    :param zarr_dir: Path to zarr files
    :param id_file: Id of file to be used from msr_files
    :param key: Key of dictionary
    :param time_diff_ref: Time difference allowed between end of experiment and end of beads measurement.
    Bead recording is shorter than experiment
    :param exclude_beads: List of BEAD IDs, e.g. ['R10', 'R01'] to exclude
    :return:
    """
    mfx = MfxData(msr_files[key][id_file], outdir_main=out_dir[key],
                          zarr_dir_main=zarr_dir[key])
    print(mfx.msrfile_path)
    mfx.MAX_TDIFF_REF = time_diff_ref
    mfx.zarr_import()
    mfx.set_valid_ref()
    if len(mfx.valid_ref_beads) == 0:
        print('No valid reference beads ')
        return
    if exclude_beads is not None:
        mfx.valid_ref_beads = [x for x in mfx.valid_ref_beads if x not in exclude_beads]
    regis = mfx.get_ref_transform()
    mfx.show_ref_transform(translate=regis[mfx.TRANS], rotate=None, save=False, show=True)
