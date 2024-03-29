import pandas as pd
import numpy as np
from sklearn import cluster
from .logformatter.customformatter import CustomFormatter
import os
from sklearn import mixture
from . import mfxcolnames as col
from evtk import hl
import logging
import sys
# Run clustering on data remove, make invalid if below a certain level separate or merge cluster


class ProcessLocalizations:
    CLS_METHOD_GMM = 'gmm'
    CLS_METHOD_BAYES_GMM = 'bayes_gmm'
    CLS_METHOD_DBSCAN = 'dbscan'

    GMM_MAX_COMPONENTS = 3          # Max components for GMM. 1 - is one fluorophore, 2 - fluorophores alternate, 3 - 2 fluorophores alternate + also simultaneously
    GMM_COVARIANCE_TYPE = 'diag'
    GMM_N_INIT = 3                  # GMM number of initializations
    STD_QUANTILE = 0.75           # STD quantile for further processing TID and split those

    def __init__(self, file_path):
        # file_path is a numpy, pickled file from  mfxdata.export_numpy()
        self.file_path = file_path

        # Some instance variables that may change depending on sample, ...
        self.MIN_LOCALIZATIONS = 3
        self.FACT_SPLIT_LOCALIZATIONS = 2     # Perform analysis and split eventually localizations with more than  self.FACT_SPLIT_LOCATIONS*self.MIN_LOCALIZATIONS
        self.DBCLUSTER_SIZE = 3               # minimal size of cluster for localizations ~ at least MIN_LOCALIZATIONS
        self.DBCLUSTER_EPS_TRACK = 2e-8       # 1e-8 = 10 nm, eps for clustering localization of a single track
        self.DBCLUSTER_EPS_MEAS = 2e-8        # dbscan eps for clustering localizations in each measurement
        self.DBCLUSTER_EPS_ALL = 2e-8         # dbscan eps for clustering localizations in all measurements
        self.DBCLUSTER_EPS_MERGED_ALL = 2e-8  # dbscan eps for clustering merged localizations in each measurement
        self.DBCLUSTER_EPS_MERGED_MEAS = 2e-8 # dbscan eps for clustering merged localizations in all measurements

        self.loc = np.load(file_path, allow_pickle=True).item()
        if 'ref_error' in self.loc.keys():
            self.ref_error = self.loc['ref_error']
            del self.loc['ref_error']

        # add clustering columns initialize with no clusters
        for label in self.loc:
            self.loc[label][col.CLS_TRACK] = np.ones(self.loc[label][col.VLD].shape)*(0) # Default is a cluster = whole track
            self.loc[label][col.CLS_MEAS] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.CLS_ALL] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.CLS_MERGED_MEAS] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.CLS_MERGED_ALL] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.TID2] = np.ones(self.loc[label][col.VLD].shape)*(-1)
        self.set_tid2()
        self.translate_to_origin()
        self.logger = logging.getLogger('processlocalizations')

        self.set_logger()

    def translate_to_origin(self, coord = col.LTR):
        # This is useful so that the lowest point is always at 0,0,0
        col_to_translate = [col.LTR, col.LNC, col.LOC, col.LRE]
        for col_name in col_to_translate:
            pos_concat = np.concatenate([self.loc[d][col_name] for d in self.loc])
            origin = np.min(pos_concat, axis = 0)
            for label in self.loc:
                self.loc[label][col_name] = self.loc[label][col_name] - origin

    def set_logger(self):
        log_fmt = "%(asctime)s [%(levelname)s] **** %(funcName)s ****\n%(message)s"
        self.logger.handlers = []
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler(self.file_path_no_ext() + ".log", mode = 'w')
        fh.setLevel(logging.INFO)
        fh.setFormatter(CustomFormatter(log_fmt))
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(CustomFormatter(log_fmt))
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)
        self.logger.propagate = False

    def log_parameters(self):
        log_msg = ""
        for k, v in self.__dict__.items():
            if k == 'loc':
                continue
            log_msg += '%s %s\n' % (k, v)
        self.logger.info(log_msg)

    def file_path_no_ext(self):
        return os.path.splitext(self.file_path)[0]

    def trim_min_localizations(self):
        """
        Trim all tracks that have not enough localizations. Set VLD to False
        """
        log_msg = ''
        for label in self.loc:
            remove_tid = 0
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                if counts[i] < self.MIN_LOCALIZATIONS:
                    remove_tid += 1
                    self.loc[label][col.VLD][range(idx, idx+counts[i])] = False
            # apply to whole data set and remove invalid localizations
            vld = self.loc[label][col.VLD]
            for col_name in self.loc[label]:
                self.loc[label][col_name] = self.loc[label][col_name][vld]
            log_msg += '%s Removed %d/%d tracks with less than %d localizations\n' % (label, remove_tid, len(uid), self.MIN_LOCALIZATIONS)
        self.logger.info(log_msg)
        # Renumerate to have continuous entries
        self.set_tid2()


    def get_var(self):
        var_all = []
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                var = np.var(self.loc[label][col.LTR][tid_idxs], axis=0)
                var_all.append(var)
        return var_all

    def get_overall_std(self, quantiles):
        var_all = self.get_var()
        std_all = np.sqrt(np.mean(var_all, axis=1))
        return {'mean': np.mean(std_all), 'quantiles': np.quantile(std_all, quantiles)}
    
    def get_split_events(self):
        keys = list(self.loc.keys())
        tid_splits = dict(zip(keys, [list()]*len(keys)))
        for label in self.loc:
            tids = []
            uid = np.unique(self.loc[label][col.TID])
            for idx in uid:
                tid_idxs = (self.loc[label][col.TID] == idx)
                uid_cls_track = np.unique(self.loc[label][col.CLS_TRACK][tid_idxs])
                if len(uid_cls_track) > 1:
                   tids.append([idx, sum(tid_idxs)])
            tid_splits[label] = tids
        return tid_splits

    def cluster_tid(self, method, coord=col.LTR, weight_prior = 1/3):
        # Cluster each track individually using heuristic of MIN_SPLIT_LOCALIZATIONS and
        # Large cluster bigger than upper certain value
        # Check if efo is not large than 2*median
        # split and clean up tracks
        # Update tid2 index
        std_all = self.get_overall_std([self.STD_QUANTILE])
        log_msg = ''
        for label in self.loc:
            efo_upper_limit = np.quantile(self.loc[label][col.EFO], [0.5])*2 # Twice higher than median, more than one fluorophore
            processed = 0
            uid, counts = np.unique(self.loc[label][col.TID], return_counts=True)
            # compute some stats, will be used to choose whether to process the track or not
            for i, uid_idx in enumerate(uid):
                tid_idxs = (self.loc[label][col.TID] == uid_idx)
                sd = np.sqrt(np.mean(np.var(self.loc[label][coord][tid_idxs], axis=0)))

                if counts[i] < self.FACT_SPLIT_LOCALIZATIONS*self.MIN_LOCALIZATIONS or sd < std_all['quantiles'][0]:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = 0
                    continue

                processed += 1
                if method == self.CLS_METHOD_DBSCAN:
                    predict = self.dbscan(self.loc[label][coord][tid_idxs],
                                          eps=self.DBCLUSTER_EPS_TRACK)
                    u_cls = np.unique(predict)
                    if len(u_cls) > 1:
                        if (len(u_cls) > 2) or (len(u_cls) == 2 and -1 not in u_cls):
                            predict = self.filter_cls_efo(predict, efo=self.loc[label][col.EFO][tid_idxs],
                                                          efo_upper_limit=efo_upper_limit)
                    self.loc[label][col.CLS_TRACK][tid_idxs] = predict
                    self.set_tid2()
                if method == self.CLS_METHOD_GMM:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.gmm(self.loc[label][coord][tid_idxs])
                    predict = self.filter_cls_efo(predict, efo=self.loc[label][col.EFO][tid_idxs],
                                                  efo_upper_limit=efo_upper_limit)
                    self.loc[label][col.CLS_TRACK][tid_idxs] = predict
                    self.set_tid2()
                if method == self.CLS_METHOD_BAYES_GMM:
                    predict = self.bayes_gmm(self.loc[label][coord][tid_idxs],
                                             std_limit=std_all['quantiles'][0],
                                             weight_prior=weight_prior)
                    predict = self.filter_cls_efo(predict, efo=self.loc[label][col.EFO][tid_idxs],
                                                  efo_upper_limit=efo_upper_limit)
                    self.loc[label][col.CLS_TRACK][tid_idxs] = predict
                    self.set_tid2()
                uid2, idxs2, counts2 = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            log_msg += '%s MIN_SPLIT_LOCALIZATION: %d, sd_limit: %.2f nm\n'  % (label, self.FACT_SPLIT_LOCALIZATIONS * self.MIN_LOCALIZATIONS, std_all['quantiles'][0] * 1e9)
            log_msg += 'Processed TID: %d / %d, Total tracks TID2: %d\n' % (processed, len(uid), len(uid2))
        self.logger.info(log_msg)

    def set_tid2(self):
        # Set index of track accounting for splitting within a track
        tid2 = 0
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for tid1 in uid:
                tid_idxs = (self.loc[label][col.TID] == tid1)
                uid_cls, idxs_cls, counts_cls = np.unique(self.loc[label][col.CLS_TRACK][tid_idxs],
                                                          return_index=True, return_counts=True)
                if len(uid_cls) == 1:
                    self.loc[label][col.TID2][tid_idxs] = tid2
                    tid2 += 1
                    continue
                else:
                    cls = self.loc[label][col.CLS_TRACK][tid_idxs]
                    tid_assign = cls.copy()
                    for i_cls in range(0, 10):
                        tid_assign[cls == i_cls] = tid2
                        tid2 = tid2 + 1
                    self.loc[label][col.TID2][tid_idxs] = tid_assign

    def cluster_meas(self, method):
        # Cluster each measurement
        # loop through labels

        for label in self.loc:
            if method == self.CLS_METHOD_DBSCAN:
                # only cluster the valid traces according to col.TID2
                idx_tid2 = self.loc[label][col.TID2] >= 0
                self.loc[label][col.CLS_MEAS][idx_tid2] = self.dbscan(self.loc[label][col.LTR][idx_tid2],
                                                                      eps=self.DBCLUSTER_EPS_MEAS)
            else:
                raise ValueError("Available methods for clustering a measurement:" + self.CLS_METHOD_DBSCAN)
        self.inter_cluster_cleanup(column_cls=col.CLS_MEAS)
        log_msg = ''
        for label in self.loc:
            uid = np.unique(self.loc[label][col.TID2])
            uid_cls = np.unique(self.loc[label][col.CLS_MEAS])
            uid_cls = uid_cls[uid_cls >= 0]
            uid = uid[uid >= 0]
            log_msg += '%s, Method %s, eps: %.2f nm\n' % (label, method, self.DBCLUSTER_EPS_MEAS*1e9)
            log_msg += 'Total tracks TID2: %d, total clusters meas: %d\n' % (len(uid), len(uid_cls))
        self.logger.info(log_msg)

    def cluster_all(self, method):
        # Concatenate the localizations
        if len(self.loc) == 1:
            self.logger.info('*******cluster_all*******\nOnly one wash. Nothing to do')
            return

        loc_all = np.concatenate([self.loc[label][col.LTR][self.loc[label][col.TID2] >= 0] for label in self.loc])
        if method == self.CLS_METHOD_DBSCAN:
            cls_idx = self.dbscan(loc_all, eps=self.DBCLUSTER_EPS_ALL)
        else:
            raise ValueError("Available methods for clustering all data:" + self.CLS_METHOD_DBSCAN)
        # Reassign per label
        istart = 0
        for label in self.loc:
            iend = istart + len(self.loc[label][col.CLS_ALL][self.loc[label][col.TID2] >= 0])
            self.loc[label][col.CLS_ALL][self.loc[label][col.TID2] >= 0] = cls_idx[range(istart, iend)]
            istart = iend

        self.inter_cluster_cleanup(column_cls=col.CLS_ALL)
        log_msg = ''
        for label in self.loc:
            uid_cls = np.unique(self.loc[label][col.CLS_ALL])
            uid_tid = np.unique(self.loc[label][col.TID2])
            uid_cls = uid_cls[uid_cls >= 0]
            uid_tid = uid_tid[uid_tid >= 0]
            log_msg += '%s, Method: %s, eps: %.2f nm\nTotal tracks TID2: %d, total cluster all: %d\n' \
                       % (label, method, self.DBCLUSTER_EPS_ALL*1e9, len(uid_tid), len(uid_cls))
        self.logger.info(log_msg)

    def cluster_all_intersect(self, column):
        # TODO: Generalize for more than 2 washes

        keys = list(self.loc.keys())
        if len(keys) == 1:
            self.logger.info("Only one wash")
            return
        if len(keys) > 2:
            self.logger.info("More than two washes. Intersect not yet implemented")
            return
        uid_P = [np.unique(self.loc[keys[0]][column]), np.unique(self.loc[keys[1]][column])]
        uid_P = [ui[ui > 0] for ui in uid_P]


        # these are the intersected joint clusters
        intersect = np.intersect1d(uid_P[0], uid_P[1])

        # count TID2 per each value intersect cluster
        tid2_intersect = []
        tid2_total = []
        for key in keys:
            idx_intersect = [cls in intersect for cls in self.loc[key][column]]
            tid2_intersect.append(np.unique(self.loc[key][col.TID2][idx_intersect]))
            tmp = np.unique(self.loc[key][col.TID2])
            tid2_total.append(tmp[tmp >= 0])
        log_msg = ''
        for idx in range(0, 2):
            k_idx1 = 0
            k_idx2 = 1
            if idx == 1:
                k_idx1 = 1
                k_idx2 = 0
            log_msg += '%s cluster with %s cls %d/%d = %.2f' % (keys[k_idx1], keys[k_idx2], len(intersect),
                                                                len(uid_P[idx]), len(intersect)/len(uid_P[idx]))
            log_msg += 'TID2 %d/%d = %.2f\n' % (len(tid2_intersect[idx]), len(tid2_total[idx]), len(tid2_intersect[idx])/len(tid2_total[idx]))
        self.logger.info(log_msg)

    def inter_cluster_cleanup(self, column_cls):
        # Loop through the data set and perform a majority voting
        # Each track can only belong to one cluster
        # Care for TID2 == -1 values

        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)

            for i, idx in enumerate(idxs):
                if uid[i] == -1:
                    continue
                tid_idxs = (self.loc[label][col.TID2] == uid[i])
                u_cls, idx_cls, counts_cls = np.unique(self.loc[label][column_cls][tid_idxs],
                                                       return_index=True, return_counts=True)
                if len(u_cls) == 1:
                    continue
                # majority vote, this could also be no clustering
                id_max = np.argsort(counts_cls)
                self.loc[label][column_cls][tid_idxs] = u_cls[id_max[-1]]

        # Assign TID2 with cluster = -1, to a unique Cluser_ID or All clusters.
        # find maximal cluster_ID
        cluster_id = np.max(np.concatenate([self.loc[label][column_cls] for label in self.loc])) + 1
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = (self.loc[label][col.TID2] == uid[i])
                if uid[i] == -1:
                    continue
                if self.loc[label][column_cls][tid_idxs][0] == -1:
                    self.loc[label][column_cls][tid_idxs] = cluster_id
                    cluster_id += 1

    def dbscan(self, x, eps, min_samples=None):

        if min_samples is None:
            min_samples = self.DBCLUSTER_SIZE
        alg = cluster.DBSCAN(eps=eps, min_samples=min_samples)
        alg.fit(x)
        pr = alg.labels_
        return pr

    def gmm(self, x):
        # run the clustering and find best BIC
        # GMM 1 component diag covariances_ gives the var in each direction
        # GMM 1 component spherical covoriances_ is the mean(var) in all directions
        lowest_bic = np.infty
        bic = []
        x = x/1e-9 # Convert to nm to have reasonable variances and not change parametrs of algortithm
        for n in range(1, self.GMM_MAX_COMPONENTS + 1):
            alg = mixture.GaussianMixture(n_components=n,
                                          covariance_type=self.GMM_COVARIANCE_TYPE,
                                          n_init=self.GMM_N_INIT)
            alg.fit(x)
            bic.append(alg.bic(x))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_alg = alg
        return best_alg.predict(x)

    def bayes_gmm(self, x, std_limit=0, weight_prior=1/3):
        x = x*1e9 # rescale in nm to simplify the fit
        std_limit = std_limit*1e9
        gmm = mixture.BayesianGaussianMixture(n_components=self.GMM_MAX_COMPONENTS,
                                              covariance_type=self.GMM_COVARIANCE_TYPE,
                                              n_init=self.GMM_N_INIT,
                                              weight_concentration_prior=weight_prior)
        ff = gmm.fit(x)
        pr = gmm.predict(x)

        if len(np.unique(pr)) == 1 or std_limit == 0:
            return pr

        u_cls, idx_cls, counts_cls = np.unique(pr, return_index=True, return_counts=True)
        # find the entry
        id_max = np.argsort(counts_cls)
        d = np.linalg.norm(gmm.means_[1] - gmm.means_[0])  # distance between centers
        if d < std_limit:
            # Merge if clouds are too close
            pr[:] = 0
        else:

            if counts_cls[id_max[0]] < self.MIN_LOCALIZATIONS:
                # Remove if cloud does not have enough points
                pr[pr == u_cls[id_max[0]]] = -1
        return pr

    def filter_cls_efo(self, pr, efo, efo_upper_limit):
        u_cls, idx_cls, counts_cls = np.unique(pr, return_index=True, return_counts=True)
        for id_cls in u_cls:
            if id_cls == -1:
                continue
            if sum(efo[pr == id_cls] > efo_upper_limit) >= round(len(efo[pr == id_cls])/2):
                pr[pr == id_cls] = -1
        return pr

    def summary_per_tid2(self):
        summary = {}
        for label in self.loc:
            summary[label] = {col.TID2: [], col.TIM: [], col.LOC: [], col.LNC: [], col.LTR: [], col.STD: [],
                              col.SE: [], col.NL: [], col.STD_XYZ: [], col.SE_XYZ: [], col.CLS_MEAS: [],
                              col.CLS_ALL: [], col.CLS_MERGED_MEAS: [], col.CLS_MERGED_ALL: []}
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                # do not add non-valid localizations where the single track failed
                if uid[i] == -1:
                    continue
                tid_idxs = (self.loc[label][col.TID2] == uid[i])
                summary[label][col.TID2].append(uid[i])
                nl = len(self.loc[label][col.LTR][tid_idxs])
                var_xyz = np.var(self.loc[label][col.LTR][tid_idxs], axis=0)
                std_xyz = np.sqrt(var_xyz)
                se_xyz = std_xyz/np.sqrt(nl)
                std = np.sqrt(np.mean(var_xyz))
                se = std/np.sqrt(nl)
                summary[label][col.TIM].append(np.mean(self.loc[label][col.TIM][tid_idxs]).tolist())
                summary[label][col.LOC].append(np.mean(self.loc[label][col.LOC][tid_idxs], axis=0).tolist())
                summary[label][col.LNC].append(np.mean(self.loc[label][col.LNC][tid_idxs], axis=0).tolist())
                summary[label][col.LTR].append(np.mean(self.loc[label][col.LTR][tid_idxs], axis=0).tolist())
                summary[label][col.STD].append(std)
                summary[label][col.SE].append(se)
                summary[label][col.NL].append(nl)
                summary[label][col.STD_XYZ].append(std_xyz)
                summary[label][col.SE_XYZ].append(se_xyz)
                summary[label][col.CLS_MEAS].append(self.loc[label][col.CLS_MEAS][tid_idxs][0])
                summary[label][col.CLS_ALL].append(self.loc[label][col.CLS_ALL][tid_idxs][0])
                summary[label][col.CLS_MERGED_MEAS].append(self.loc[label][col.CLS_MERGED_MEAS][tid_idxs][0])
                summary[label][col.CLS_MERGED_ALL].append(self.loc[label][col.CLS_MERGED_ALL][tid_idxs][0])




                # cluster single wash measurement at the level of merged localizations
        for label in summary:
            summary[label][col.CLS_MERGED_MEAS] = self.dbscan(summary[label][col.LTR],
                                                              eps=self.DBCLUSTER_EPS_MERGED_MEAS, min_samples=1)


        # cluster combined wash at the level of merged localizations
        if len(summary.keys()) > 1:
            loc_all = np.concatenate([summary[label][col.LTR] for label in summary])
            cls_idx = self.dbscan(loc_all, eps=self.DBCLUSTER_EPS_MERGED_ALL, min_samples=1)
            istart = 0
            for label in summary:
                iend = istart + len(summary[label][col.CLS_ALL])
                summary[label][col.CLS_MERGED_ALL] = cls_idx[range(istart, iend)]
                istart = iend
        return summary

    @staticmethod
    def filter_summary_loc(in_dict, filter_cluster):
        out_dict = in_dict.copy()
        for label in in_dict:
            uid, counts = np.unique(in_dict[label][filter_cluster['name']], return_counts=True)
            uid_valid = uid[counts >= filter_cluster['size']]
            idx_valid = np.where(np.in1d(in_dict[label][filter_cluster['name']], uid_valid))[0]
            for col_name in in_dict[label]:
                out_dict[label][col_name] = np.array(in_dict[label][col_name])[idx_valid.astype(int)]
        return out_dict

    def export_csv(self, in_dict, file_path):
        print('export csv table')
        keys = list(in_dict.keys())
        for i, key in enumerate(keys):
            df = pd.DataFrame.from_dict(in_dict[key], orient='columns')
            # expand coordinates entries
            for colname in [col.LTR, col.LOC, col.LNC, col.STD_XYZ, col.SE_XYZ]:
                df[[colname + '_x', colname + '_y', colname + '_z']] = pd.DataFrame(
                    df[colname].tolist(), index=df.index)
            df['wash'] = np.repeat(key, len(df[colname].tolist()))
            if i == 0:
                df.to_csv(file_path + '.csv', mode='w', index=False)
            else:
                df.to_csv(file_path + '.csv', mode='a', index=False, header=False)

    def export_raw_vtu(self, in_dict, coord, file_path, unit=col.UNIT_M):
        # TODO: Filter by cluster size. Typically clusters of size <2 are singletons
        print('export paraview file')
        # export one single file
        export_column = [col.TID2, col.TIM]
        # export one file per wash

        if unit == col.UNIT_M:
            factor_unit = 1
        if unit == col.UNIT_UM:
            factor_unit = 1e6
        if unit == col.UNIT_NM:
            factor_unit = 1e9

        label_idx = 0
        for label in in_dict:
            label_idx += 1
            # concatenate positions
            pos_concat = np.array(in_dict[label][coord])*factor_unit
            x = np.ascontiguousarray(pos_concat[:, 0], dtype=np.float64)
            y = np.ascontiguousarray(pos_concat[:, 1], dtype=np.float64)
            z = np.ascontiguousarray(pos_concat[:, 2], dtype=np.float64)
            pview_dict = {}
            for column in export_column:
                if column == col.SE:
                    pview_dict[column] = np.array(in_dict[label][column])*factor_unit
                else:
                    pview_dict[column] = np.array(in_dict[label][column])

            pview_dict['wash'] = np.repeat(label_idx, len(in_dict[label][coord]))
            hl.pointsToVTK(file_path + '_' + coord + '_' + label, x, y, z,
                           data=pview_dict)
        if len(in_dict) > 1:
            # concatenate all washes
            # export to vtk format for view in paraview, ideally also clustering etc.
            pos = np.concatenate([in_dict[d][coord] for d in in_dict])*factor_unit
            x = np.ascontiguousarray(pos[:, 0], dtype=np.float64)
            y = np.ascontiguousarray(pos[:, 1], dtype=np.float64)
            z = np.ascontiguousarray(pos[:, 2], dtype=np.float64)
            pview_dict = {}
            for column in export_column:
                if column == col.SE:
                    pview_dict[column] = np.array([in_dict[d][column] for d in in_dict])*factor_unit
                else:
                    pview_dict[column] = np.concatenate([in_dict[d][column] for d in in_dict])
            keys = list(in_dict.keys())
            pview_dict['wash'] = np.repeat(range(1, len(keys)+1), [len(in_dict[k][coord]) for k in keys])
            hl.pointsToVTK(file_path + '_' + coord, x, y, z, data=pview_dict)

    def export_vtu(self, in_dict, coord, file_path, unit=col.UNIT_M):
        # TODO: Filter by cluster size. Typically clusters of size <2 are singletons
        print('export paraview file')
        # export one single file
        export_column = [col.TID2, col.TIM, col.CLS_ALL, col.CLS_MEAS,
                         col.CLS_MERGED_MEAS, col.CLS_MERGED_ALL, col.SE, col.NL]
        # export one file per wash

        if unit == col.UNIT_M:
            factor_unit = 1
        if unit == col.UNIT_UM:
            factor_unit = 1e6
        if unit == col.UNIT_NM:
            factor_unit = 1e9

        label_idx = 0
        for label in in_dict:
            label_idx += 1
            # concatenate positions
            pos_concat = np.array(in_dict[label][coord])*factor_unit
            x = np.ascontiguousarray(pos_concat[:, 0], dtype=np.float64)
            y = np.ascontiguousarray(pos_concat[:, 1], dtype=np.float64)
            z = np.ascontiguousarray(pos_concat[:, 2], dtype=np.float64)
            pview_dict = {}
            for column in export_column:
                if column == col.SE:
                    pview_dict[column] = np.array(in_dict[label][column])*factor_unit
                else:
                    pview_dict[column] = np.array(in_dict[label][column])

            pview_dict['wash'] = np.repeat(label_idx, len(in_dict[label][coord]))
            pview_dict['se_norm'] = pview_dict['se']/np.max(pview_dict['se'])
            pview_dict['nl_norm'] = pview_dict['nl']/np.max(pview_dict['nl'])
            hl.pointsToVTK(file_path + '_' + coord + '_' + label, x, y, z,
                           data=pview_dict)
        if len(in_dict) > 1:
            # concatenate all washes
            # export to vtk format for view in paraview, ideally also clustering etc.
            pos = np.concatenate([in_dict[d][coord] for d in in_dict])*factor_unit
            x = np.ascontiguousarray(pos[:, 0], dtype=np.float64)
            y = np.ascontiguousarray(pos[:, 1], dtype=np.float64)
            z = np.ascontiguousarray(pos[:, 2], dtype=np.float64)
            pview_dict = {}
            for column in export_column:
                if column == col.SE:
                    pview_dict[column] = np.concatenate([in_dict[d][column] for d in in_dict])*factor_unit
                else:
                    pview_dict[column] = np.concatenate([in_dict[d][column] for d in in_dict])
            keys = list(in_dict.keys())
            pview_dict['wash'] = np.repeat(range(1, len(keys)+1), [len(in_dict[k][coord]) for k in keys])
            pview_dict['se_norm'] = pview_dict['se']/np.max(pview_dict['se'])
            pview_dict['nl_norm'] = pview_dict['nl']/np.max(pview_dict['nl'])
            hl.pointsToVTK(file_path + '_' + coord, x, y, z, data=pview_dict)







