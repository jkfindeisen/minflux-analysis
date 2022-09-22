# -----
col.CLS_MERGED_MEAS = 'cls_merged_meas'
col.CLS_MEAS = 'cls_meas'
col.CLS_MERGED_ALL = 'cls_merged_all'
col.CLS_ALL = 'cls_all'


cluster_counts <- function(afile, cls_type){
  data <-read.csv(afile)
  hist_counts <- data %>% count(data[cls_type])
  return(hist_counts)
}

cluster_counts_set <- function(files, cls_type){
  hist_counts_all <- list()
  for (afile in files) {
    hist_counts_all <- rbind(hist_counts_all, cluster_counts(afile, cls_type = cls_type))
    
  }
  return(hist_counts_all)
}
  



cluster_stats <- function(afile, cls_type){
  
  data <-read.csv(afile)
  cls_stats <- cluster_counts(afile, cls_type)
  
  cls_stats$longest_radius_nm <- 0
  cls_stats$PCA1_nm <- 0
  cls_stats$PCA2_nm <- 0
  cls_stats$PCA3_nm <- 0
  cls_stats$PCA_vol_um3 <- 0
  
  
  for (irow in seq(1, nrow(cls_stats))) {
    # to compute object you need at least 3 points > 3
    
    if(cls_stats$n[irow] < 2){
      next
    }
    clust_coord = data[data[cls_type]==cls_stats[irow, cls_type], c('ltr_x', 'ltr_y', 'ltr_z') ]
    PCA <- prcomp(clust_coord)
    dd <- dist(clust_coord)
    max_radius <- max(dd)/2
    if (cls_stats$n[irow] == 2) 
      cls_stats[irow, c('PCA1_nm', 'PCA2_nm')] <- PCA$sdev*1e9
    else {
      cls_stats[irow, c('PCA1_nm', 'PCA2_nm', 'PCA3_nm')] <- PCA$sdev*1e9
      cls_stats[irow, 'PCA_vol_um3'] <- 4/3*pi*PCA$sdev[1]*PCA$sdev[2]*PCA$sdev[3]*1e6*1e6*1e6
    }
    cls_stats[irow, 'longest_radius_nm'] <- max_radius*1e9
  }
  return(cls_stats)
}

cluster_stats_set <- function(files, cls_type){
  cls_stats_all <- list()
  for (afile in files) {
    cls_stats_all <- rbind(cls_stats_all, cluster_stats(afile, cls_type = cls_type))
    
  }
  return(cls_stats_all)
}

