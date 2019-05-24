library(venneuler)
vd <- venneuler(c("USPSTF"=2.4,"Risk-based"=0.9,"Life gained"=0.5,"USPSTF&Risk-based"=0.4,"USPSTF&Life gained"=0.7,"Risk-based&Life gained"=2.2,"USPSTF&Risk-based&Life gained"=4.8))
vd$labels <- c("","","")
vd$colors[1] <- .15
vd$colors[2] <- .9
vd$colors[3] <- .55

#colors: .1 (orange), .2(tan), .3(green), .4(green), .5(blue-green), 
#.6 (light-blue), .7 (purplish-blue), .8 (light-purple), .9 (pepto-bismal), 1 (pink)

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas5.jpeg",width=8,height=8,units='in',res=900)
plot(vd,main="Coverage of strategies selecting ever-smokers for lung cancer screening")
text(x=.815,y=0.74,"USPSTF-based")
text(x=.45,y=0.815,"Life gained based")
text(x=.43,y=0.18,"Risk-based")

text(x=.5,y=.54,"4.8m screened (8.2%)",cex=.65)
text(x=.5,y=.525,"35.2k deaths prevented (43.2%)",cex=.65)
text(x=.5,y=.51,"396k years gained (33.4%)",cex=.65)

text(x=.585,y=.27,"0.4m screened (0.6%)",cex=.65)
text(x=.585,y=.255,"1.3k deaths prevented (1.6%)",cex=.65)
text(x=.585,y=.24,"13k years gained (1.1%)",cex=.65)

text(x=.62,y=.75,"0.7m screened (1.1%)",cex=.65)
text(x=.62,y=.735,"1.5k deaths prevented (1.8%)",cex=.65)
text(x=.62,y=.72,"31k years gained (2.6%)",cex=.65)

text(x=.245,y=.575,"2.2m screened (3.7%)",cex=.65)
text(x=.245,y=.56,"10.6k deaths prevented (13.0%)",cex=.65)
text(x=.245,y=.545,"124k years gained (10.5%)",cex=.65)

text(x=.29,y=.76,"0.5m screened (0.9%)",cex=.65)
text(x=.29,y=.745,"1.2k deaths prevented (1.4%)",cex=.65)
text(x=.29,y=.73,"23k years gained (2.0%)",cex=.65)

text(x=.25,y=.3,"0.9m screened (1.6%)",cex=.65)
text(x=.25,y=.285,"3.1k deaths prevented (3.9%)",cex=.65)
text(x=.25,y=.27,"30k years gained (2.6%)",cex=.65)

text(x=.795,y=.455,"2.4m screened (4.1%)",cex=.65)
text(x=.795,y=.44,"3.4k deaths prevented (4.2%)",cex=.65)
text(x=.795,y=.425,"64k years gained (5.4%)",cex=.65)

text(x=.8,y=.195,"47m unscreened (80%)",cex=.65)
text(x=.8,y=.18,"25.2k preventable deaths (30.9%)",cex=.65)
text(x=.8,y=.165,"502k years lost (42.4%)",cex=.65)
dev.off()