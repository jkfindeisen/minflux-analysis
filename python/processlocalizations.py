
import numpy as np
from sklearn import cluster
import os
from sklearn import mixture
import mfxcolnames as col

# Run clustering on data remove, make invalid if below a certain level separate or merge cluster

class ProcessLocalizations:
    CLS_METHOD_GMM = 'gmm'
    CLS_METHOD_BAYES_GMM = 'bayes_gmm'
    CLS_METHOD_DBSCAN = 'dbscan'

    GMM_MAX_COMPONENTS = 2
    GMM_COVARIANCE_TYPE = 'diag'
    GMM_N_INIT = 3                  # GMM number of initializations
    MIN_LOCALIZATIONS = 4           # Minimum number of localizations pre track
    MIN_WEIGHT_GMM = 0.4            # Baysian mixture model to allow for 2 Gaussians
    MIN_SPLIT_LOCALIZATIONS = 2*MIN_LOCALIZATIONS    # Perform analysis and split eventually localizations with more than 8
    DBCLUSTER_EPS_TRACK = 1e-8    # 1e-8 = 10 nm, eps for clustering localization of a single track
    DBCLUSTER_EPS_MEAS = 3e-8     # eps for clustering localization of all localization in one measurement
    DBCLUSTER_EPS_ALL = 3e-8      # eps for clustering localization of all localization in all measurement
    DBCLUSTER_SIZE = 4

    def __init__(self, file_path):
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
        remove_tid = 0
        total_tid = 0
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            total_tid += len(uid)
            for i, idx in enumerate(idxs):
                if counts[i] < self.MIN_LOCALIZATIONS:
                    remove_tid += 1
                    self.loc[label][col.VLD][range(idx, idx+counts[i])] = False
            # apply to whole data set and remove invalid localizations
            vld = self.loc[label][col.VLD]
            for col_name in self.loc[label]:
                self.loc[label][col_name] = self.loc[label][col_name][vld]
        # Renumerate to have continuous entries
        self.set_tid2()
        print("Removed %d out of %d tracks with less than %d localisations" % (remove_tid, total_tid, self.MIN_LOCALIZATIONS))

    def get_var(self):
        var_all = []
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                var = np.var(self.loc[label][col.LTR][tid_idxs], axis=0)
                var_all.append(var)
        return var_all

    def get_overall_std(self):
        var_all = self.get_var()
        std_all = np.sqrt(np.mean(var_all, axis=1))
        return {'mean': np.mean(std_all), 'quantiles': np.quantile(std_all, [0.25, 0.5, 0.75])}

    def cluster_tid(self, method):
        # Cluster each track individually using
        # loop through labels
        std_all = self.get_overall_std()
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            # compute some stats, will be used to choose whether to process the track or not
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                sd = np.sqrt(np.mean(np.var(self.loc[label][col.LTR][tid_idxs], axis=0)))
                if counts[i] < self.MIN_SPLIT_LOCALIZATIONS or sd < std_all['quantiles'][2]:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = 0
                    continue
                if method == self.CLS_METHOD_DBSCAN:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.dbscan(self.loc[label][col.LTR][tid_idxs],
                                                                           eps=self.DBCLUSTER_EPS_TRACK)
                    self.set_tid2()
                if method == self.CLS_METHOD_GMM:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.gmm(self.loc[label][col.LTR][tid_idxs])
                    self.set_tid2()
                if method == self.CLS_METHOD_BAYES_GMM:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.bayes_gmm(self.loc[label][col.LTR][tid_idxs],
                                                                              std_limit=std_all['quantiles'][2])
                    self.set_tid2()

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

    def cluster_all(self, method):
        # Concatenate the localizations

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

    def inter_cluster_cleanup(self, column_cls):
        # Loop through the data set and perform a majority voting
        # Each track can only belong to one cluster
        # Remove -1 entries
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID2], return_index=True, return_counts=True)

            for i, idx in enumerate(idxs):
                if uid[i] == -1:
                    continue
                tid_idxs = range(idx, idx+counts[i])
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
        print(cluster_id)

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

    def summary_per_tid2(self):
        summary = {}
        for label in self.loc:
            summary[label] = {col.TID2: [], col.LTR:[], col.STD: [], col.SE: [], col.NL: [],
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
                mean_pos = np.mean(self.loc[label][col.LTR][tid_idxs], axis=0)

                summary[label][col.LTR].append(mean_pos.tolist())
                summary[label][col.STD].append(std)
                summary[label][col.SE].append(se)
                summary[label][col.NL].append(nl)
                summary[label][col.STD_XYZ].append(std_xyz)
                summary[label][col.SE_XYZ].append(se_xyz)
                if self.loc[label][col.CLS_MEAS][tid_idxs[0]] < 0:
                    print(self.loc[label][col.CLS_MEAS][tid_idxs])
                summary[label][col.CLS_MEAS].append(self.loc[label][col.CLS_MEAS][tid_idxs[0]])
                summary[label][col.CLS_ALL].append(self.loc[label][col.CLS_ALL][tid_idxs[0]])
        return summary


    def bayes_gmm(self, x, std_limit=0):
        x = x/1e-9 # rescale in nm to simplify the fit
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
                # Remove if cloud does not have anough points
                pr[pr == u_cls[id_max[0]]] = -1
        return pr
    # def get_gmmpar(self, gmm, indata):
    #     if gmm.covariance_type == "full":
    #         covariances = gmm.covariances_
    #     elif gmm.covariance_type == "tied":
    #         covariances = gmm.covariances_[:2, :2]
    #     elif gmm.covariance_type == "diag":
    #         sd = np.sqrt(gmm.covariances_)
    #         m = gmm.means_
    #         v = None
    #         return [m, v]
    #     elif gmm.covariance_type == "spherical":
    #         covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_
    #     v, w = np.linalg.eigh(covariances)
    #     u = w[0] / np.linalg.norm(w[0])
    #     angle = np.arctan2(u[1], u[0])
    #     angle = 180 * angle / np.pi  # convert to degrees
    #     v = 2.0 * np.sqrt(2.0) * np.sqrt(v) # confidence interval ellipse

if __name__ == "__main__":
    os.environ["OMP_NUM_THREADS"] = '1'
    pl = ProcessLocalizations(
        "C:/Users/apoliti/Desktop/mfluxtest/analysis/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.npy")
    print(pl.DBCLUSTER_EPS_TRACK)
    pl.trim_min_localizations()
    pl.cluster_tid(method=pl.CLS_METHOD_BAYES_GMM)
    pl.cluster_meas(method=pl.CLS_METHOD_DBSCAN)
    pl.cluster_all(method=pl.CLS_METHOD_DBSCAN)
    summary = pl.summary_per_tid2()


    std_all = pl.get_overall_std()
    print(std_all)



