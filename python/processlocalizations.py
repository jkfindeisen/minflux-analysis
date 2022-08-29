from mfxdata import MfxData
import numpy as np


# Run clustering on data remove, make invalid if below a certain level separate or merge cluster

class ProcessLocalizations:
    MIN_LOCALIZATIONS = 3
    SPLIT = True
    MERGE = True

    def __init__(self, file_path):
        self.loc = np.load(file_path, allow_pickle=True).item()

    def min_localizations(self):
        removeTid = 0
        totalTid = 0
        for label in self.loc:
            for tid in np.unique(self.loc[label]['tid']):
                totalTid += 1
                idxs = np.where(self.loc[label]['tid'] == tid)[0]
                if len(idxs) < self.MIN_LOCALIZATIONS:
                    removeTid += 1
                    self.loc[label]['vld'][idxs] = False
        print("Removed %d out of %d tracks with less than %d localisations" % (removeTid, totalTid, self.MIN_LOCALIZATIONS))

    def merge_tracks(self):
        
if __name__ == "__main__":
    pl = ProcessLocalizations(
        "C:/Users/apoliti/Desktop/mfluxtest/analysis/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b/220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b.npy")
    print(pl.loc)
    pl.min_localizations()

