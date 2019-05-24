for (i in 1:5){
  nhis <- data.frame()
  for (j in 2013:2015) {
    dirname <- paste0("/home/cheunglc/lyg/nhis2013_15", i, j)
    load(file=paste0(dirname,"/y_mod.RData"))
    nhis <- rbind(nhis,y)
  }
  save(nhis,file=paste0("/home/cheunglc/lyg/nhis2013_15/nhis_imputed_mod", i,".RData"))
}


