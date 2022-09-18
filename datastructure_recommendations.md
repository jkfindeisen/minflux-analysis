# Expected data structure (!project specific!)

These recommendations apply to acquisitions linked by repeated images on the same structure (DNA-paint washes).
To simplify the processing part of the code expects a specific format of the data. 

The data structure may differ from what has been acquired originally on the microscope. Where typically all ROIS
are within a msr file. 

From the original file one should create in *Imspector* additional files 
that fulfill the expected format.

## Directory structure

It is preferred that the raw data (*msr* files)  is kept in a different folder than the analysis results.
One directory per experiment, e.g. named by the date YYYYMMDD, or very clear naming of the files that
allow to identify the files that belong to the same experiment.

## msr file name

Very consistent name, that contains a unique identifier for the experiment (e.g. the date), some
information on the sample and the washes.

## msr file content

* It is preferred that the raw data (* \*.msr *) is kept in a different folder than the analysis results
* One single *msr* file per ROI containing all subsequent washes.
* The Minflux data is named `P1`, `P2`, etc. corresponding to the order of the washes. In *Imspector* this is done in `Workspace > Minflux data > Right-mouse-click on data > Rename`
* `README.txt`file and/or excel `README.xls` file that contains relevan

# msr file name
Please follow these rules for the data folder:

    Strict separation of raw data from analysis results.
    Remove all mat, npy files or anything that is not raw data.
    In the raw data directory ONE msr file per ROI. No merge of ROIs in one file. See Naming of msr. Remove all files that do not fulfill this and eventually duplicate the data.
    The single msr file contains all washes. In Impspector > Worksapce > Minflux  the wahses are (re)named P1, P2, ... . The name correspond to the order in which the washes have been performed.
    If needed create a README.txt file or excel file with important comments. Like wash order, protein, laser settings, pinhole size.
    If there has been a problem, like lost beads, or reassignment of the beads due to a restart please specify it here.  

Naming-scheme of msr file:

The key is consistency. I would like to browse through the data and process it in an automatic way. This is not possible for the moment. Please rename all files accordingly!

These are example names I found in your folders

    220309_VGlut_paint_2nM_3DMINFLUX_16%_PH0,6_P1_03.msr
    220510_ATG_paint_2nM_3DMINFLUX_16p_PH0,6_P2_03b.msr
    220601_Syp_P1_ATG9_P2_ROI_01.msr
    220825_ATG9_paint_2nM_3DMINFLUX_16p_PH06_P2_Roi_1.msr
    220804_VGlut1_NB_paint_2nM_3DMINFLUX_16p_PH06_P1_ROI1_first30min.msr
    220510_Syp_1.msr

Here is my suggestion:

Short version:


    YYMMDD_{name-P1}_{name-P2}_ROI{roiID}.msr

example:

    220510_ATG9_Syp_ROI03.msr

The {name-P1}, ... must be consistent. ATG9 or ATG? VGlut or VGlut1?

The {roiID} matches your original, during acquisition saved record. This version scales better with more washes. The rest of the information should be in a README.txt file + labbook,...

Long version (or anything you think is important)

    YYMMDD_{minflux-method}_{name-P1}_{technique-P1}_{name-P2}_{technique-P2}_{laser-settings?}p-PH{pinhole-size}_ROI{roiID}.msr

example

    220510_3DMINFLUX_ATG9_paint_2nM_Syp_paint_2nM_16p_PH06_ROI03.msr

This format is too long and redundant. Some imaging settings can also be stored in an excel table that belongs in the directory of the data.



