for (i in 1:5){
  nhis <- data.frame()
  for (j in 2006:2009) {
       dirname <- paste0("/home/cheunglc/lyg/nhis2006_09", i, j)
       load(file=paste0(dirname,"/y_mod.RData"))
       nhis <- rbind(nhis,y)
  }
  save(nhis,file=paste0("/home/cheunglc/lyg/nhis2006_09_mod", i,".RData"))
}
