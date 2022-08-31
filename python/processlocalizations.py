
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
    GMM_N_INIT = 3
    MIN_LOCALIZATIONS = 3           # Minimum number of localizations pre track
    MIN_WEIGHT_GMM = 0.4            # Baysian mixture model to allow for 2 Gaussians
    MIN_SPLIT_LOCALIZATIONS = 10    # Perform analysis and split eventually localizations with more than
    SPLIT = True
    MERGE = True
    DBCLUSTER_EPS_TRACK = 1e-8    # 1e-8 = 10 nm, eps for clustering localization of a single track
    DBCLUSTER_EPS_MEAS = 3e-8     # eps for clustering localization of all localization in one measurement
    DBCLUSTER_EPS_ALL = 3e-8      # eps for clustering localization of all localization in all measurement
    DBCLUSTER_SIZE = 3

    def __init__(self, file_path):
        self.loc = np.load(file_path, allow_pickle=True).item()
        # add clustering columns initialize with no clusters
        for label in self.loc:
            self.loc[label][col.CLS_TRACK] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.CLS_MEAS] = np.ones(self.loc[label][col.VLD].shape)*(-1)
            self.loc[label][col.CLS_ALL] = np.ones(self.loc[label][col.VLD].shape)*(-1)

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
            # apply to whole data set
            vld = self.loc[label][col.VLD]
            for col_name in self.loc[label]:
                self.loc[label][col_name] = self.loc[label][col_name][vld]
        print("Removed %d out of %d tracks with less than %d localisations" % (remove_tid, total_tid, self.MIN_LOCALIZATIONS))

    def get_stats(self):
        std_all = []
        for label in self.loc:
            vld = self.loc[label]['vld']
            unique_tid = np.unique(self.loc[label]['tid'][vld])
            for tid in unique_tid:
                idxs = np.where(self.loc[label]['tid'] == tid)[0]
                std = np.std(self.loc[label]['ltr'][idxs], axis=0)
                #print(np.linalg.norm(std))
                std_all.append(std)
        return std_all

    def cluster_tid(self, method):
        # Cluster each track individually using
        # loop through labels
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                if method == self.CLS_METHOD_DBSCAN:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.dbscan(self.loc[label][col.LTR][tid_idxs], eps=self.DBCLUSTER_EPS_TRACK)
                if method == self.CLS_METHOD_GMM:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.gmm(self.loc[label][col.LTR][tid_idxs])
                if method == self.CLS_METHOD_BAYES_GMM:
                    self.loc[label][col.CLS_TRACK][tid_idxs] = self.bayes_gmm(self.loc[label][col.LTR][tid_idxs])

    def cluster_meas(self, method):
        # Cluster each measurement
        # loop through labels
        for label in self.loc:
            if method == self.CLS_METHOD_DBSCAN:
                self.loc[label][col.CLS_MEAS] = self.dbscan(self.loc[label][col.LTR], eps=self.DBCLUSTER_EPS_MEAS)
            else:
                raise ValueError("Available methods for clustering a measurement:" + self.CLS_METHOD_DBSCAN)

    def cluster_all(self, method):
        # Concatenate the localizations

        loc_all = np.concatenate([self.loc[label][col.LTR] for label in self.loc])
        if method == self.CLS_METHOD_DBSCAN:
            cls_idx = self.dbscan(loc_all, eps=self.DBCLUSTER_EPS_ALL)
        else:
            raise ValueError("Available methods for clustering all data:" + self.CLS_METHOD_DBSCAN)
        # Reassign per label
        istart = 0
        for label in self.loc:
            iend = istart + len(self.loc[label][col.CLS_ALL])
            self.loc[label][col.CLS_ALL] = cls_idx[range(istart, iend)]
            istart = iend

    def cluster_cleanup(self, column_cls):
        # Loop through the data set and perform a majority voting to which tracks
        # Kepp the -1 entries
        for label in self.loc:
            uid, idxs, counts = np.unique(self.loc[label][col.TID], return_index=True, return_counts=True)
            for i, idx in enumerate(idxs):
                tid_idxs = range(idx, idx+counts[i])
                u_cls, idx_cls, counts_cls = np.unique(self.loc[label][column_cls][tid_idxs], return_index=True, return_counts=True)
                print(u_cls)
                if len(u_cls) == 1:
                    continue
                if len(u_cls) == 2 and any(c == -1 for c in u_cls):
                    continue
                # majority vote
                non_clust = np.where(self.loc[label][column_cls][tid_idxs] == -1)
                id_max = np.argsort(counts_cls)
                if u_cls[id_max[-1]] != -1:
                    self.loc[label][column_cls][tid_idxs] = u_cls[id_max[-1]]
                else:
                    self.loc[label][column_cls][tid_idxs] = u_cls[id_max[-2]]
                # Reassign non_clustered values
                for i_non_clust in non_clust[0].tolist():
                    self.loc[label][column_cls][tid_idxs[i_non_clust]] = -1

    def bayes_gmm(self, x):
        x = x/1e-9
        gmm = mixture.BayesianGaussianMixture(n_components=self.GMM_MAX_COMPONENTS,
                                              covariance_type=self.GMM_COVARIANCE_TYPE,
                                              n_init=self.GMM_N_INIT)
        ff = gmm.fit(x)
        return gmm.predict(x)

    def dbscan(self, x, eps):
        alg = cluster.DBSCAN(eps=eps, min_samples=self.DBCLUSTER_SIZE)
        alg.fit(x)
        return alg.labels_

    def gmm(self, x):
        # run the clustering and find best BIC
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

    def get_gmmpar(self, gmm, indata):
        if gmm.covariance_type == "full":
            covariances = gmm.covariances_
        elif gmm.covariance_type == "tied":
            covariances = gmm.covariances_[:2, :2]
        elif gmm.covariance_type == "diag":
            sd = np.sqrt(gmm.covariances_)
            m = gmm.means_
            v = None
            return [m, v]
        elif gmm.covariance_type == "spherical":
            covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_
        v, w = np.linalg.eigh(covariances)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)

if __name__ == "__main__":
    os.environ["OMP_NUM_THREADS"] = '1'
    pl = ProcessLocalizations(
        "C:/Users/apoliti/Desktop/mfluxtest/analysis/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.npy")
    print(pl.DBCLUSTER_EPS_TRACK)
    pl.trim_min_localizations()
    #pl.cluster_tid(method=pl.CLS_METHOD_DBSCAN)
    #pl.cluster_meas(method=pl.CLS_METHOD_DBSCAN)
    pl.cluster_all(method=pl.CLS_METHOD_DBSCAN)
    pl.cluster_cleanup(col.CLS_ALL)
    pl.cluster_tid(method=pl.CLS_METHOD_BAYES_GMM)



