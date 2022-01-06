calc_svds_fun<-function(
  cxy,
  left_field,
  right_field,
  loc_lf,
  loc_rf,
  n
){
  # Notification:
  # DO follow the configuration guidance before using the function
  # Do put this function into the same directory with gpu_mca_analysis.R

  library(data.table)
  library(bigstatsr)
  library(rARPACK)
  
  cxy = as.matrix(cxy)
  system.time(d1 <- svds(cxy,k = n))
  
  ds = diag(d1$d)
  u = d1$u
  v = d1$v
  
  fwrite(ds,'ds.csv') # divergence
  fwrite(u,'u_ecx.csv') # spatial modes of left field
  fwrite(v,'v_ecy.csv') # spatial modes of right field
  
  # input x and y 
  
  calc_ais_fun <- function(i){
    ui = u[,i]
    tmp_ai = t(ui) %*% left_field
    return(tmp_ai)
  }
  
  calc_bis_fun <- function(i){
    vi = v[,i]
    tmp_bi = t(vi) %*% right_field
    return(tmp_bi)
  }
  
  
  i = 1:ncol(u)
  a_mat = do.call(rbind,lapply(i,calc_ais_fun)) # rows to models
  b_mat = do.call(rbind,lapply(i,calc_bis_fun))
  
  fwrite(a_mat,'ajs_full_mat.csv') # temporal mode of left field
  fwrite(b_mat,'bjs_full_mat.csv') # temporal mode of right field
  
  calc_covs_between_ai_bi<-function(i){
    tmp_ai = a_mat[i,]
    tmp_bi = b_mat[i,]
    
    tmp_covs = cov(tmp_ai,tmp_bi)
    return(tmp_covs)
  }
  
  i = 1:nrow(a_mat)
  
  # explained covariance for each mode
  covs_bet_ai_bi = do.call('c',lapply(i,calc_covs_between_ai_bi))
  
  modu = u[,1:4] # 1st to 4th modes
  modv = v[,1:4]
  
 
  colnames(loc_lf) = c('long','lat')
  colnames(loc_rf) = c('long','lat')
  
  
  colnames(modu) = paste0('Model_left',1:4)
  colnames(modv) = paste0('Model_right',1:4)
  
  modu = data.frame(loc_lf,modu)
  modv = data.frame(loc_rf,modv)
  
  library(ggplot2)
  
  modu = reshape2::melt(modu,c('long','lat'))
  modv = reshape2::melt(modv,c('long','lat'))
  
  minu = min(modu$value)
  maxu = max(modu$value)
  stepu = (maxu - minu)/9
  
  minv = min(modv$value)
  maxv = max(modv$value)
  stepv = (maxv - minv)/9
  
  modu$levels = cut(modu$value,
                    breaks = seq(minu,maxu,stepu))
  modv$levels = cut(modv$value,
                    breaks = seq(minv,maxv,stepv))
  
  modu$type = rep(paste0('Model ',1:4),each = nrow(loc_lf))
  modv$type = rep(paste0('Model ',1:4),each = nrow(loc_rf))
  
  
  fwrite(modu,'modu.csv') # 1st to 4th spatial modes of left field
  fwrite(modv,'modv.csv') # 1st to 4th spatial modes of right field
  
  
  
  p1 = ggplot()+
    geom_tile(data = modu,
              aes(x = long,y = lat,fill = levels))+
    #scale_fill_distiller(palette = 'Spectral')+
    scale_fill_brewer(palette = 'Spectral')+
    facet_wrap(~variable,nrow = 4)+
    theme_bw()+
    theme(legend.position = 'bottom')
  
  p2 = ggplot()+
    geom_tile(data = modv,
              aes(x = long,y = lat,fill = levels))+
    #scale_fill_distiller(palette = 'Spectral')+
    scale_fill_brewer(palette = 'Spectral')+
    facet_wrap(~variable,nrow = 4)+
    theme_bw()+
    theme(legend.position = 'bottom')
  
  library(cowplot)
  p12 = plot_grid(
    p1,p2,nrow = 1,
    rel_widths = c(1,1),
    rel_heights = c(1,1)
  )
  
  png('test_mca.png',
      height = 25,
      width = 27,
      units = 'cm',
      res = 800)
  print(p12)
  dev.off()
  
  
  
}