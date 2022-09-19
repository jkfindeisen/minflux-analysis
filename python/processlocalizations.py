
import numpy as np
from sklearn import cluster
import os
from sklearn import mixture
import mfxcolnames as col
from evtk import hl

# Run clustering on data remove, make invalid if below a certain level separate or merge cluster

class ProcessLocalizations:
    CLS_METHOD_GMM = 'gmm'
    CLS_METHOD_BAYES_GMM = 'bayes_gmm'
    CLS_METHOD_DBSCAN = 'dbscan'

    GMM_MAX_COMPONENTS = 3          # Max components for GMM. 1 - is one fluorophore, 2 - fluorophores alternate, 3 - 2 fluorophores alternate + also simultaneously
    GMM_COVARIANCE_TYPE = 'diag'
    GMM_N_INIT = 3                  # GMM number of initializations
    MIN_LOCALIZATIONS = 4           # Minimum number of localizations pre track
    MIN_SPLIT_LOCALIZATIONS = 2*MIN_LOCALIZATIONS    # Perform analysis and split eventually localizations with more than 8
    DBCLUSTER_EPS_TRACK = 2e-8    # 1e-8 = 10 nm, eps for clustering localization of a single track
    DBCLUSTER_EPS_MEAS = 2e-8     # eps for clustering localization of all localization in one measurement
    DBCLUSTER_EPS_ALL = 2e-8      # eps for clustering localization of all localization in all measurement
    DBCLUSTER_SIZE = 4
    STD_QUANTILE = 0.75           # STD quantile for further processing TID and split those

    def __init__(self, file_path):
        # file_path is a pickled file
        self.loc = np.load(file_path, allow_pickle=True).item()
        # add clustering columns initialize with no clusters
        for label in self.loc:
            self.loc[label][col.CLS_TRACK] = np.ones(self.loc[label][col.VLD].shape)*(0) # Default is a cluster = whole track
            self.loc[label][col.CLS_MEAS] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.CLS_ALL] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.TID2] = np.ones(self.loc[label][col.VLD].shape)*(-1)
        self.set_tid2()

    def trim_min_localizations(self):
        """
        Trim all tracks that have not enough localizations. Set VLD to False
        """

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
            print("*******time_min_localizations*******\n%s Removed %d/%d tracks with less than %d localisations"
                  % (label, remove_tid, len(uid), self.MIN_LOCALIZATIONS))
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

    def cluster_tid(self, method):
        # Cluster each track individually using heuristic of MIN_SPLIT_LOCALIZATIONS and
        # Large cluster bigger than upper certain value
        # Check if efo is not large than 2*median
        # split and clean up tracks
        # Update tid2 index
        std_all = self.get_overall_std([self.STD_QUANTILE])

        for label in self.loc:
            efo_upper_limit = np.quantile(self.loc[label][col.EFO], [0.5])*2 # Twice higher than median, more than one fluorophore

            processed = 0
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            # compute some stats, will be used to choose whether to process the track or not
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                sd = np.sqrt(np.mean(np.var(self.loc[label][col.LTR][tid_idxs], axis=0)))

                if counts[i] < self.MIN_SPLIT_LOCALIZATIONS or sd < std_all['quantiles'][0]:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = 0
                    continue

                processed += 1
                if method == self.CLS_METHOD_DBSCAN:
                    predict = self.dbscan(self.loc[label][col.LTR][tid_idxs],
                                          eps=self.DBCLUSTER_EPS_TRACK)
                    u_cls = np.unique(predict)
                    if len(u_cls) > 1:
                        if (len(u_cls) > 2) or (len(u_cls) == 2 and -1 not in u_cls):
                            predict = self.filter_cls_efo(predict, efo=self.loc[label][col.EFO][tid_idxs],
                                                          efo_upper_limit=efo_upper_limit)
                    self.loc[label][col.CLS_TRACK][tid_idxs] = predict
                    self.set_tid2()
                if method == self.CLS_METHOD_GMM:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.gmm(self.loc[label][col.LTR][tid_idxs])
                    predict = self.filter_cls_efo(predict, efo=self.loc[label][col.EFO][tid_idxs],
                                                  efo_upper_limit=efo_upper_limit)
                    self.loc[label][col.CLS_TRACK][tid_idxs] = predict
                    self.set_tid2()
                if method == self.CLS_METHOD_BAYES_GMM:
                    predict = self.bayes_gmm(self.loc[label][col.LTR][tid_idxs],
                                             std_limit=std_all['quantiles'][0])
                    predict = self.filter_cls_efo(predict, efo=self.loc[label][col.EFO][tid_idxs],
                                                  efo_upper_limit=efo_upper_limit)
                    self.loc[label][col.CLS_TRACK][tid_idxs] = predict
                    self.set_tid2()
                uid2, idxs2, counts2 = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            print('*******cluster_tid*******\n%s MIN_SPLIT_LOCALIZATION: %d, sd_limit: %.2f nm\nProcessed TID: %d / %d, Total tracks TID2: %d'
                  % (label, self.MIN_SPLIT_LOCALIZATIONS, std_all['quantiles'][0]*1e9, processed, len(uid), len(uid2)))

    def set_tid2(self):
        # Set index of track accounting for splitting within a track
        tid = 0
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = list(range(idx, idx+counts[i]))
                uid_cls, idxs_cls, counts_cls = np.unique(self.loc[label][col.CLS_TRACK][tid_idxs],
                                                          return_index=True, return_counts=True)
                if len(uid_cls) == 1:
                    self.loc[label][col.TID2][tid_idxs] = tid
                    tid += 1
                    continue
                else:
                    cls = self.loc[label][col.CLS_TRACK][tid_idxs]
                    tid_assign = cls.copy()
                    for i_cls in range(0, 10):
                        tid_assign[cls == i_cls] = tid
                        tid = tid + 1
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
        for label in self.loc:
            uid = np.unique(self.loc[label][col.TID2])
            uid_cls = np.unique(self.loc[label][col.CLS_MEAS])
            uid_cls = uid_cls[uid_cls >= 0]
            uid = uid[uid >=0 ]
            print('*******cluster_meas*******\n%s, Method %s, eps: %.2f nm\nTotal tracks TID2: %d, total clusters meas: %d'
                  % (label, method, self.DBCLUSTER_EPS_MEAS*1e9, len(uid), len(uid_cls)))

    def cluster_all(self, method):
        # Concatenate the localizations
        if len(self.loc) == 1:
            print('*******cluster_all*******\nOnly one wash. Nothing to do')
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
        for label in self.loc:
            uid_cls = np.unique(self.loc[label][col.CLS_ALL])
            uid_tid = np.unique(self.loc[label][col.TID2])
            uid_cls = uid_cls[uid_cls >= 0]
            uid_tid = uid_tid[uid_tid >= 0]
            print('*******cluster_all*******\n%s, Method: %s, eps: %.2f nm' % (label, method, self.DBCLUSTER_EPS_ALL*1e9))
            print('Total tracks TID2: %d, total cluster all: %d' % (len(uid_tid), len(uid_cls)))

    def cluster_all_intersect(self):
        # TODO: Generalize for more than 2 washes
        keys = list(self.loc.keys())
        uid_P = [np.unique(self.loc[keys[0]][col.CLS_ALL]), np.unique(self.loc[keys[1]][col.CLS_ALL])]
        uid_P = [ui[ui > 0] for ui in uid_P]


        # these are the intersected joint clusters
        intersect = np.intersect1d(uid_P[0], uid_P[1])

        # count TID2 per each value intersect cluster
        tid2_intersect = []
        tid2_total = []
        for key in keys:
            idx_intersect = [cls in intersect for cls in self.loc[key][col.CLS_ALL]]
            tid2_intersect.append(np.unique(self.loc[key][col.TID2][idx_intersect]))
            tmp = np.unique(self.loc[key][col.TID2])
            tid2_total.append(tmp[tmp >= 0])

        print('*******cluster_all_intersect*******')
        for idx in range(0, 2):
            k_idx1 = 0
            k_idx2 = 1
            if idx == 1:
                k_idx1 = 1
                k_idx2 = 0
            print('%s cluster with %s cls %d/%d = %.2f; TID2 %d/%d = %.2f' %
                  (keys[k_idx1], keys[k_idx2], len(intersect), len(uid_P[idx]), len(intersect)/len(uid_P[idx]),
                   len(tid2_intersect[idx]), len(tid2_total[idx]), len(tid2_intersect[idx])/len(tid2_total[idx])))

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

        # Assign TID2 with cluster = -1, to a unique Cluser_ID or All clusters
        # find maximal cluster_ID
        cluster_id = np.max(np.concatenate([self.loc[label][column_cls] for label in self.loc])) + 1
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                if uid[i] == -1:
                    continue
                if self.loc[label][column_cls][idx] == -1:
                    self.loc[label][column_cls][tid_idxs] = cluster_id
                    cluster_id += 1

    def dbscan(self, x, eps):
        alg = cluster.DBSCAN(eps=eps, min_samples=self.DBCLUSTER_SIZE)
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

    def bayes_gmm(self, x, std_limit=0, efo = None, efo_upper_limit = None):
        x = x*1e9 # rescale in nm to simplify the fit
        std_limit=std_limit*1e9
        gmm = mixture.BayesianGaussianMixture(n_components=self.GMM_MAX_COMPONENTS,
                                              covariance_type=self.GMM_COVARIANCE_TYPE,
                                              n_init=self.GMM_N_INIT)
        ff = gmm.fit(x)
        pr = gmm.predict(x)

        if len(np.unique(pr)) == 1 or std_limit == 0:
            return pr

        u_cls, idx_cls, counts_cls = np.unique(pr, return_index=True, return_counts=True)
        # find the entry
        id_max = np.argsort(counts_cls)
        d = np.linalg.norm(gmm.means_[1] - gmm.means_[0]) # distance between centers
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
        # TODO: perform additional clustering from average distrubtion
        summary = {}
        for label in self.loc:
            summary[label] = {col.TID2: [], col.TIM: [], col.LTR: [], col.STD: [], col.SE: [], col.NL: [],
                              col.STD_XYZ: [], col.SE_XYZ: [],  col.CLS_MEAS: [], col.CLS_ALL: []}
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                if uid[i] == -1: # do not add non-valid localizations
                    continue
                tid_idxs = range(idx, idx+counts[i])
                summary[label][col.TID2].append(uid[i])
                nl = len(self.loc[label][col.LTR][tid_idxs])
                var_xyz = np.var(self.loc[label][col.LTR][tid_idxs], axis=0)
                std_xyz = np.sqrt(var_xyz)
                se_xyz = std_xyz/np.sqrt(nl)
                std = np.sqrt(np.mean(var_xyz))
                se = std/np.sqrt(nl)
                mean_time = np.mean(self.loc[label][col.TIM][tid_idxs])
                mean_pos = np.mean(self.loc[label][col.LTR][tid_idxs], axis=0)
                summary[label][col.TIM].append(mean_time)
                summary[label][col.LTR].append(mean_pos.tolist())
                summary[label][col.STD].append(std)
                summary[label][col.SE].append(se)
                summary[label][col.NL].append(nl)
                summary[label][col.STD_XYZ].append(std_xyz)
                summary[label][col.SE_XYZ].append(se_xyz)
                # track ID where clustering did not succeed
                if self.loc[label][col.CLS_MEAS][tid_idxs[0]] < 0:
                    print('CLS_MEAS failed for tif %d', self.loc[label][col.CLS_MEAS][tid_idxs])
                summary[label][col.CLS_MEAS].append(self.loc[label][col.CLS_MEAS][tid_idxs[0]])
                summary[label][col.CLS_ALL].append(self.loc[label][col.CLS_ALL][tid_idxs[0]])
        return summary

    def export_vtu(self, lcoord, file_path):
        out_dict = self.summary_per_tid2()
        # export to vtk format for view in paraview, ideally also clustering etc.

        pos_concat = np.concatenate([out_dict[d][lcoord] for d in out_dict])
        tid_concat = np.concatenate([out_dict[d][col.TID2] for d in out_dict])
        tim_concat = np.concatenate([out_dict[d][col.TIM] for d in out_dict])
        cls_all = np.concatenate([out_dict[d][col.CLS_ALL] for d in out_dict])
        cls_meas = np.concatenate([out_dict[d][col.CLS_MEAS] for d in out_dict])
        se = np.concatenate([out_dict[d][col.SE] for d in out_dict])
        nl = np.concatenate([out_dict[d][col.NL] for d in out_dict])

        # concatenate positions
        x = np.ascontiguousarray(pos_concat[:, 0], dtype = np.float64)
        y = np.ascontiguousarray(pos_concat[:, 1], dtype = np.float64)
        z = np.ascontiguousarray(pos_concat[:, 2], dtype = np.float64)
        #se_xyz = np.ascontiguousarray(se_xyz)

        keys = list(out_dict.keys())
        lp1_p2 = [len(out_dict[k][lcoord]) for k in keys]
        p = np.repeat(range(1, len(keys)+1), [len(out_dict[k][lcoord]) for k in keys])
        # No string output

        hl.pointsToVTK(file_path, x, y, z, data={'wash': p, col.TID2: tid_concat,
                                                   col.CLS_ALL: cls_all,
                                                   col.CLS_MEAS: cls_meas,
                                                   col.SE: se/np.min(se),
                                                   col.NL: nl, col.TIM: tim_concat
                                                 })
        for label in out_dict:
            # concatenate positions
            pos_concat = np.array(out_dict[label][lcoord])
            x = np.ascontiguousarray(pos_concat[:, 0], dtype=np.float64)
            y = np.ascontiguousarray(pos_concat[:, 1], dtype=np.float64)
            z = np.ascontiguousarray(pos_concat[:, 2], dtype=np.float64)

            hl.pointsToVTK(file_path + label, x, y, z,
                           data={col.TID2: np.array(out_dict[label][col.TID2]),
                                 col.TIM: np.array(out_dict[label][col.TIM]),
                                 col.SE: np.array(out_dict[label][col.SE]/np.min(se)),
                                 col.NL: np.array(out_dict[label][col.NL]/np.min(nl)),
                                 col.CLS_ALL: np.array(out_dict[label][col.CLS_ALL]),
                                 col.CLS_MEAS: np.array(out_dict[label][col.CLS_MEAS])
                                 })

if __name__ == "__main__":
    os.environ["OMP_NUM_THREADS"] = '1'
    pl = ProcessLocalizations(
        'Z:/siva_minflux/analysis//Multiwash/Syp_ATG9/220510_Syp_ATG9_ROI01\\220510_Syp_ATG9_ROI01.npy')
    pl.trim_min_localizations()
    pl.cluster_tid(method=pl.CLS_METHOD_BAYES_GMM)
    pl.cluster_meas(method=pl.CLS_METHOD_DBSCAN)
    pl.cluster_all(method=pl.CLS_METHOD_DBSCAN)
    pl.cluster_all_intersect()
    pl.export_vtu(col.LTR, 'Z:/siva_minflux/analysis//Multiwash/Syp_ATG9/220510_Syp_ATG9_ROI01\\220510_Syp_ATG9_ROI01_ltr')
    #summary = pl.summary_per_tid2()


    #std_all = pl.get_overall_std()
    #print(std_all)



