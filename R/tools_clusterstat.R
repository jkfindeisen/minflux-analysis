# -----
col.CLS_MERGED_MEAS = 'cls_merged_meas'
col.CLS_MEAS = 'cls_meas'
col.CLS_MERGED_ALL = 'cls_merged_all'
col.CLS_ALL = 'cls_all'
protname_synonyms <- data.frame('Syp' = 'Syp', 'SYP' = 'Syp', 'Syn' = 'Syp',
                                'ATG9' = 'ATG9', 
                                'ZnT3' = 'ZnT3', 
                                'VGLUT1' = 'VGLUT1', 
                                'Picc' = 'Picc')

cluster_counts <- function(afile, cls_type){
  data <-read.csv(afile)
  pattern <- '(^\\d+)(-||_)([[:alnum:]]+)(_P||-P||--P)(\\d)'
  
  # check the number of washes
  wash <- unique(data$wash)
  hist_counts_all <- list()
  
  for (awash in wash) {
    exp_date <- gsub(pattern = pattern, awash, replacement = '\\1')
    protname <- gsub(pattern = pattern, awash, replacement = '\\3')
    id_wash <-  gsub(pattern = pattern, awash, replacement = '\\5')
    hist_counts <- data[data$wash == awash, ] %>% count(data[data$wash == awash, cls_type])
    colnames(hist_counts) <- c(cls_type, 'n')
    hist_counts$wash <- awash
    hist_counts$protname <- protname_synonyms[protname][[1]]
    hist_counts$exp_date <- exp_date
    hist_counts$id_wash <- id_wash
    hist_counts_all <- rbind(hist_counts_all, hist_counts)
  }
  if (cls_type %in% c(col.CLS_ALL, col.CLS_MERGED_ALL)){
    
    if (length(wash) == 2) {
      if (cls_type == col.CLS_ALL) {
        hist_counts_p1p2 <- hist_counts_all %>% group_by(cls_all) %>% 
          summarise(n_all = sum(n), 
                    P1 = sum(n[id_wash==1]),
                    P2 = sum(n[id_wash==2]))
      }
      if (cls_type == col.CLS_MERGED_ALL) {
        hist_counts_p1p2 <- hist_counts_all %>% group_by(cls_merged_all) %>% 
          summarise(n_all = sum(n), 
                    P1 = sum(n[id_wash==1]),
                    P2 = sum(n[id_wash==2]) )
        
      }
      hist_counts_all <- hist_counts_p1p2
      hist_counts_all$protname <- NA
      hist_counts_all$exp_date <- exp_date
      colnames(hist_counts_all)[colnames(hist_counts_all) == 'n_all'] <- 'n'
    }
  }
  return(hist_counts_all)
}

cluster_counts_set <- function(files, cls_type){
  hist_counts_all <- list()
  for (afile in files) {
    hist_counts_all <- rbind(hist_counts_all, cluster_counts(afile, cls_type = cls_type))
    
  }
  return(hist_counts_all)
}

cluster_size <- function(data, cls_stats, cls_type){
  for (irow in seq(1, nrow(cls_stats))) {
    tryCatch({
      # to compute object you need at least 3 points > 3
      clust_coord <- data[data[cls_type]==cls_stats[irow, cls_type][[1]], c('ltr_x', 'ltr_y', 'ltr_z') ]
      cls_stats[irow, c('centroid_x', 'centroid_y', 'centroid_z')] <-
        as.list(c(mean(clust_coord$ltr_x), mean(clust_coord$ltr_y), mean(clust_coord$ltr_z)))
      if(cls_stats$n[irow] < 2){
        next
      }
      
      PCA <- prcomp(clust_coord)
      dd <- dist(clust_coord)
      max_radius <- max(dd)/2
      if (cls_stats$n[irow] == 2) 
        cls_stats[irow, c('PCA1_nm', 'PCA2_nm')] <- as.list(PCA$sdev[1:2]*1e9)
      else {
        cls_stats[irow, c('PCA1_nm', 'PCA2_nm', 'PCA3_nm')] <- as.list(PCA$sdev*1e9)
        cls_stats[irow, 'PCA_vol_um3'] <- 4/3*pi*PCA$sdev[1]*PCA$sdev[2]*PCA$sdev[3]*1e6*1e6*1e6
      }
      cls_stats[irow, 'longest_radius_nm'] <- max_radius*1e9
    }, warning=function(w) print(irow), error=function(e) print(paste0('Error ', irow, ' ', e)))
  }
  return(cls_stats)
  
}

cluster_stats <- function(afile, cls_type){
  data <- read.csv(afile)
  cls_stats <- cluster_counts(afile, cls_type)
  cls_stats$longest_radius_nm <- 0
  cls_stats$PCA1_nm <- 0
  cls_stats$PCA2_nm <- 0
  cls_stats$PCA3_nm <- 0
  cls_stats$PCA_vol_um3 <- 0
  cls_stats$centroid_x <- 0
  cls_stats$centroid_y <- 0
  cls_stats$centroid_z <- 0
  cls_stats$file <- basename(dirname(afile))
  dim_img <- apply(cbind(data$ltr_x, data$ltr_y, data$ltr_x), 2, max) - apply(cbind(data$ltr_x, data$ltr_y, data$ltr_x), 2, min)
  cls_stats$dim_img_vol_m3 <- prod(dim_img)
  
  # process each wash independently to compute controids etc
  if (cls_type %in% c(col.CLS_MEAS, col.CLS_MERGED_MEAS)) {
    for (awash in unique(cls_stats$wash)){
      datasub <- subset(data, wash == awash)
      cls_stats_sub <- subset(cls_stats, wash == awash)
      cls_stats[cls_stats$wash == awash, ] <- cluster_size(data = datasub, cls_stats = cls_stats_sub, cls_type = cls_type)
    }
  } 
  if (cls_type %in% c(col.CLS_ALL, col.CLS_MERGED_ALL)) {
    cls_stats <- cluster_size(data = data, cls_stats = cls_stats, cls_type = cls_type)
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

# display a cluster and different metrics


show_cluster <- function(afile, cls_type, threeD = TRUE) {
  
  data <-read.csv(afile)
  cls_counts <- cluster_counts(afile, cls_type)
  # find biggest cluster
  cls_counts <- subset(cls_counts, id_wash == 2 & n >=5)
  cls_counts <- cls_counts[order(cls_counts$n, decreasing = TRUE), ]
  data <- subset(data, wash == cls_counts$wash[1])
  rand_id <- sample(seq(2,length(cls_counts$n)), size = 1)
  
  big_cls <- cls_counts[rand_id, cls_type]
  clust_coord <- data[data[cls_type]==big_cls, c('ltr_x', 'ltr_y', 'ltr_z') ]
  PCA <- prcomp(clust_coord, scale = FALSE)

  dd <- dist(clust_coord)
  max_size <- max(dd)
  idx_max_elong <- which(as.matrix(dd) == max_size, arr.ind = TRUE)[1,]
  dt <- cbind(x = clust_coord$ltr_x, y = clust_coord$ltr_y, z =  clust_coord$ltr_z)
  ellipse <- rgl::ellipse3d(cov(dt), centre = PCA$center, level = 0.68)
  # Create coordinate pairs to display  
  x <- list()
  y <- list()
  z <- list()
  for (id_vec in c(1,2,3)){
    
    x <- rbind(x, cbind(as.vector(PCA$center[1] - PCA$rotation[1,id_vec]*PCA$sdev[id_vec]), as.vector(PCA$center[1] + PCA$rotation[1,id_vec]*PCA$sdev[id_vec])))
    y <- rbind(y, cbind(as.vector(PCA$center[2] - PCA$rotation[2,id_vec]*PCA$sdev[id_vec]), as.vector(PCA$center[2] + PCA$rotation[2, id_vec]*PCA$sdev[id_vec])))
    z <- rbind(z, cbind(as.vector(PCA$center[3] - PCA$rotation[3,id_vec]*PCA$sdev[id_vec]), as.vector(PCA$center[3] + PCA$rotation[3, id_vec]*PCA$sdev[id_vec])))
  }
  x <- rbind(x, cbind(clust_coord[idx_max_elong[1], 1],clust_coord[idx_max_elong[2], 1]))
  y <- rbind(y, cbind(clust_coord[idx_max_elong[1], 2],clust_coord[idx_max_elong[2], 2]))
  z <- rbind(z, cbind(clust_coord[idx_max_elong[1], 3],clust_coord[idx_max_elong[2], 3]))
  scene <- list(camera = list(eye = list(x = 1.25, y = 0, z = 1.25)), 
                xaxis = list(title='x'), 
                yaxis = list(title ='y'), 
                zaxis=list(title = 'z'), 
                aspectmode = 'data', 
                aspectratio=list(x=1, y=1, z=1))
  
  fig <- plot_ly() %>% add_trace(data=data[data[cls_type]==big_cls, ], x = ~ltr_x, y = ~ltr_y, z = ~ltr_z, #, size=~se, sizes=c(5, 50), 
                                 marker = list(symbol = 'diamond', sizemode = 'diameter'), 
                                 mode='markers', type='scatter3d',  name=paste('Cluster', big_cls), showlegend =TRUE)
  
  fig <- fig %>% add_trace(x=~x[1,] , y=~y[1,], z=~z[1,],  type="scatter3d", mode="lines", line = list(width=10), 
                             name=paste('Diameter ', round(2*PCA$sdev[1]*1e9, 1), '(nm)'), showlegend = TRUE ) 
  #fig <- fig %>% add_trace(x=~x[4,] , y=~y[4,], z=~z[4,], type="scatter3d", mode="lines", line = list(width=10), 
  #                         name=paste('Maximal extension ', round(max_size*1e9, 1), '(nm)'),  
  #                         showlegend = TRUE)
  fig <- fig %>% layout(scene=scene)
  # ----
  if (!threeD) {
    print(max_size)
    lbl = paste( round(2*PCA$sdev[1]*1e9, 1), 'nm')
    print(lbl)
  fig <- autoplot(PCA, loadings = FALSE, geom = 'point',  scale =0, shape = 'square', size = 4, colour = 'gray')
  fig <- fig + geom_segment(aes(x = -PCA$sdev[1], xend = PCA$sdev[1], 
                            y =0 , yend = 0), color = 'blue', size = 1.5)
  
  
  fig <- fig + theme_bw()
  fig <- fig + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        panel.grid = element_blank())
  fig <- fig + scale_x_continuous(limits = c(-2.5*PCA$sdev[1],2.5*PCA$sdev[1] ))
  fig <- fig + scale_y_continuous(limits = c(-2.5*PCA$sdev[1],2.5*PCA$sdev[1] ))
  #fig <- fig + annotate('text', x =0, y = 3e-9, parse = FALSE,
  #                      label = paste0('Diameter 2*sigma ',  lbl )
  #                        )
  fig
  }
  # ----
  return(fig)
}

# set.seed(6)
# set.seed(6)
# cls_counts <- cluster_counts(csv_files[1], cls_type = col.CLS_MERGED_MEAS)
# fig = show_cluster(csv_files[1], cls_type = col.CLS_MERGED_MEAS, threeD = FALSE)
# fig
ggsave(plot = fig, filename= 'Z:/siva_minflux/analysis/Multiwash/ZnT3_Syp/figures/clst_syp_merged_meas.png',
      width = 5 , height = 5, unit = 'cm')
