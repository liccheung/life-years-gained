for (i in 1:5){
  nhis <- data.frame()
  for (j in 2002:2009) {
    dirname <- paste0("/home/cheunglc/lyg/nhis2002_09", i, j)    
    load(file=paste0(dirname,"/y_mod.RData"))
    nhis <- rbind(nhis,y)
  }
  save(nhis,file=paste0("/home/cheunglc/lyg/nhis2002_09/nhis_imputed_mod", i,".RData"))
}
