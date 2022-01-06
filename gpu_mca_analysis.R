gpu_mca_analysis <- function(
  working_path = '', # working directory in your computer
  input_leftfield = '', #spatiotemporal matrix for left field
  input_rightfield = '',#spatiotemporal matrix for right field
  input_locleft = '', # longitude and latitude for left field
  input_locright = '', # longitude and latitude for right field
  step = 20000 # computation step up to computation capability of computer GPU
               # especially up to the memory of the graphic 
){
  # Notification:
  # DO follow the configuration guidance before using the function


  library(rARPACK)
  library(bigstatsr)
  library(data.table)
  library(gpuR)
  library(gputools)
  
  setwd(working_path)
  
  #import left and right fields
  right_field = as.data.frame(fread(input_rightfield))
  loc_rf = as.data.frame(fread(input_locright))
  
  left_field = as.data.frame(fread(input_leftfield))
  loc_lf = as.data.frame(fread(input_locleft))
  
  # transfer left field to m*n spatiotemporal matrix
  # m refers to number of the spatial cells for left field
  # n refers to temporal length n
  left_field = t(left_field) 
  
  # transfer right field to q*n spatiotemporal matrix
  # q refers to number of the spatial cells for right field
  # n refers to temporal length n
  right_field = t(right_field) 
  
  nrow_left = nrow(left_field)
  nrow_right = nrow(right_field)
  n = ncol(left_field)
  
  if(nrow_left > nrow_right){
    sid = seq(1,nrow_left,step)
    eid = seq(step,nrow_left,step)
    eid = c(eid,nrow_left)
    
    cc
    
  }else{
    sid = seq(1,nrow_right,step)
    eid = seq(step,nrow_right,step)
    eid = c(eid,nrow_right)
    
    cxy = list()
    for(i in 1:length(sid)){
      tmpsid = sid[i]
      tmpeid = eid[i]
      
      right_field1 = as.matrix(right_field[tmpsid:tmpeid,])
      left_field1 = as.matrix(left_field[1:nrow_left,])
      
      system.time(
        tmpcxy <- 1/n * gpuMatMult(left_field1,t(right_field1))
      )
      cxy[[i]] = tmpcxy
      print(i)
    }
    cxy = do.call(cbind,cxy)
    
  }
  
  output_cxy= 'cxy.csv'
  fwrite(cxy,output_cxy)
  
  
  source("calc_svds_fun.R")
  # note: do place both functions into the same working directory
  # working_path here
  calc_svds_fun(
    cxy,
    left_field,
    right_field,
    loc_lf,
    loc_rf,
    n
  )
  
}