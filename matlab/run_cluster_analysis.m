[coloc_clusters, size_clusters] = minflux_colocalization_analysis({'Z:\siva_minflux\Two color\ZnT3-Syp\09-03-22\MSR files\MAT files\03\220309-151621_3Dminflux_VGlut_Paint_2nM_16%640_PH06_P1.mat',  ...
    'Z:\siva_minflux\Two color\ZnT3-Syp\09-03-22\MSR files\MAT files\03\220309-143513_3Dminflux_VGlut_Paint_2nM_16%640_PH06_P2.mat'});

%%

n1 = coloc_clusters(:,1);
n2 = coloc_clusters(:,2);

tbl = table(n1, n2, size_clusters);

writetable(tbl, "Z:\siva_minflux\example_export.xlsx")