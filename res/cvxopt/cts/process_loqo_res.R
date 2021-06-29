setwd("C:/Users/Uikos/Dropbox/vanderbei - kevin research/CS - code/res")

el1_values = c(2,50,100,125,150)
el1sparse_values = c(2,10,20,30,40,50,60,70,80,90,100,125,150)

el1_time = c(748.2,712.5,749.8,767.6,859.6,782.9,783.0,783.9,806.7,782.3)
el1sparse_time = c(47.0,42.3,44.4,44.6,46.9,49.1,44.4,49.1,46.6,48.9,46.5,48.9,46.7,49.2,48.8,49.1,48.9,48.8,51.0,51.1,48.8,51.1,53.2,52.9,55.4,57.6)


el1sparse_resnames = rep(NA,2*length(el1sparse_values))
el1sparse_resfilenames = rep(NA,2*length(el1sparse_values))
el1sparse_l1error = rep(NA,2*length(el1sparse_values))
el1sparse_maxerror = rep(NA,2*length(el1sparse_values))
el1sparse_sparsity = rep(NA,2*length(el1sparse_values)) #to check our tests are valid

el1_resnames = rep(NA,2*length(el1_values))
el1_resfilenames = rep(NA,2*length(el1_values))
el1_l1error = rep(NA,2*length(el1_values))
el1_maxerror = rep(NA,2*length(el1_values))
el1_sparsity = rep(NA,2*length(el1_values)) #to check our tests are valid

for(i in 1:length(el1sparse_values)){
  for(j in 0:1){
    el1sparse_resnames[2*i+j-1] = paste("el1sparse_res",el1sparse_values[i],"_",j,sep="")
    el1sparse_resfilenames[2*i+j-1] = paste("loqo_res_",el1sparse_values[i],"_",j,"_el1sparse.csv",sep="")
    assign(el1sparse_resnames[2*i+j-1],read.csv(file=el1sparse_resfilenames[2*i+j-1],sep="",header=FALSE))
    
    x = get(el1sparse_resnames[2*i+j-1])$V3
    y = get(el1sparse_resnames[2*i+j-1])$V4
    
    len = dim(el1sparse_res2_0)[1]
    x = x[-len]
    y = y[-len]
    
    el1sparse_l1error[2*i+j-1] = sum(abs(x-y))
    el1sparse_maxerror[2*i+j-1] = max(abs(x-y))
    el1sparse_sparsity[2*i+j-1] = sum(x)
  }
}

datatable = rbind(el1sparse_sparsity,el1sparse_time,el1sparse_maxerror,el1sparse_l1error)
datatable = t(datatable)
write.table(datatable,file="res_loqo_el1sparse.csv",sep=",",col.names=FALSE,row.names=FALSE)


for(i in 1:length(el1_values)){
  for(j in 0:1){
    el1_resnames[2*i+j-1] = paste("el1_res",el1_values[i],"_",j,sep="")
    el1_resfilenames[2*i+j-1] = paste("loqo_res_",el1_values[i],"_",j,"_el1.csv",sep="")
    assign(el1_resnames[2*i+j-1],read.csv(file=el1_resfilenames[2*i+j-1],sep="",header=FALSE))
    
    x = get(el1_resnames[2*i+j-1])$V2
    y = get(el1_resnames[2*i+j-1])$V3
    
    len = dim(el1_res2_0)[1]
    x = x[-len]
    y = y[-len]
    
    el1_l1error[2*i+j-1] = sum(abs(x-y))
    el1_maxerror[2*i+j-1] = max(abs(x-y))
    el1_sparsity[2*i+j-1] = sum(x)
  }
}

datatable = rbind(el1_sparsity,el1_time,el1_maxerror,el1_l1error)
datatable = t(datatable)
write.table(datatable,file="res_loqo_el1.csv",sep=",",col.names=FALSE,row.names=FALSE)









