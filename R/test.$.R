load(file="mortality.model.v2.RData")

for (i in 1:5){
  load(file=paste0("/home/cheunglc/lyg/nhis2002_09/nhis_imputed_", i,".RData"))
  for (j in 2002:2009) {
    eval(parse(text=paste0("y <- subset(nhis,year==",j,")")))
    
    dirname <- paste0("/home/cheunglc/lyg/nhis2002_09", i, j)
    com <- paste("mkdir", dirname)
    system(com)
    save(y,file=paste0(dirname,"/y.RData"))
    
    #create R file
    write(paste0("load('",dirname,"/y.RData')"),paste0(dirname,"/lyg_",i,j,".R"))
    file.append(paste0(dirname,"/lyg_",i,j,".R"),"/home/cheunglc/lyg/lyg.R")
    com <- paste0("save(y,file='",dirname,"/y_mod.RData')")
    cat(com,file=paste0(dirname,"/lyg_",i,j,".R"), append=TRUE)
    
    #create swarm file
    command <- paste0("R --vanilla <",dirname,"/lyg_",i,j,".R >",dirname,"/lyg_",i,j,".out")
    write(command,paste0(dirname,"/swarm"))
    
    #submit swarm file
    com<-paste0("swarm -g 1 -f ",dirname,"/swarm"," --module R/3.3.2_gcc-6.2.0", " --time 3-00:00:00")
    system(com)
    
  }
}

