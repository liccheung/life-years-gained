library(venneuler)
vd <- venneuler(c("USPSTF"=2.4,"Risk-based"=0.9,"Life gained"=0.5,
                  "USPSTF&Risk-based"=0.4,"USPSTF&Life gained"=0.7,"Risk-based&Life gained"=2.2,
                  "USPSTF&Risk-based&Life gained"=4.8))
vd$labels <- c("","","")
vd$colors[1] <- .15
vd$colors[2] <- .9
vd$colors[3] <- .55


vd2 <- venneuler(c("USPSTF"=3.4,"Risk-based"=3.1,"Life gained"=1.2,
                   "USPSTF&Risk-based"=1.3,"USPSTF&Life gained"=1.5,"Risk-based&Life gained"=10.6,
                   "USPSTF&Risk-based&Life gained"=35.2))
vd2$labels <- c("","","")
vd2$colors[1] <- .15
vd2$colors[2] <- .9
vd2$colors[3] <- .55

vd3 <- venneuler(c("USPSTF"=64,"Risk-based"=30,"Life gained"=23,
                   "USPSTF&Risk-based"=13,"USPSTF&Life gained"=31,"Risk-based&Life gained"=124,
                   "USPSTF&Risk-based&Life gained"=396))
vd3$labels <- c("","","")
vd3$colors[1] <- .15
vd3$colors[2] <- .9
vd3$colors[3] <- .55

#colors: .1 (orange), .2(tan), .3(green), .4(green), .5(blue-green), 
#.6 (light-blue), .7 (purplish-blue), .8 (light-purple), .9 (pepto-bismal), 1 (pink)

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas10.jpeg",width=11.5,height=5,units='in',res=900)
par(mfrow=c(1,3), mar=c(6,1,5,2))
plot(vd,main="Number of ever-smokers")
text(x=.815,y=0.75,"USPSTF")
text(x=.44,y=0.83,"Life gained based")
text(x=.43,y=0.17,"Risk-based")

text(x=.5,y=.54,"4.8m\n(8.2%)",cex=.75)
text(x=.585,y=.27,"0.4m\n(0.6%)",cex=.75)
text(x=.60,y=.735,"0.7m\n(1.1%)",cex=.75)
text(x=.26,y=.575,"2.2m\n(3.7%)",cex=.75)
text(x=.37,y=.76,"0.5m\n(0.9%)",cex=.75)
text(x=.27,y=.3,"0.9m\n(1.6%)",cex=.75)
text(x=.795,y=.455,"2.4m\n(4.1%)",cex=.75)
text(x=.8,y=.195,"47m unscreened\n(80%)",cex=.75)

plot(vd2,main="Number of preventable deaths")
text(x=.81,y=0.74,"USPSTF")
text(x=.45,y=0.8275,"Life gained based")
text(x=.35,y=0.1675,"Risk-based")

text(x=.5,y=.525,"35.2k\n(43.2%)",cex=.75)
arrows(x0=.72, y0=.255, x1 = .685, y1 = .28, length = 0.025, angle = 30,
       code = 2)
text(x=.75,y=.255,"1.3k\n(1.6%)",cex=.75)
arrows(x0=.65, y0=.79, x1 = .62, y1 = .765, length = 0.025, angle = 30,
       code = 2)
text(x=.65,y=.82,"1.5k\n(1.8%)",cex=.75)
text(x=.23,y=.56,"10.6k\n(13.0%)",cex=.75)
arrows(x0=.23, y0=.77, x1 = .275, y1 = .74, length = 0.025, angle = 30,
       code = 2)
text(x=.2,y=.8,"1.2k\n(1.4%)",cex=.75)
arrows(x0=.2, y0=.225, x1 = .27, y1 = .27, length = 0.025, angle = 30,
       code = 2)
text(x=.2,y=.2,"3.1k\n(3.9%)",cex=.75)
text(x=.81,y=.46,"3.4k\n(4.2%)",cex=.75)
text(x=.75,y=.18,"25.2k unprevented deaths\n(30.9%)",cex=.75)

plot(vd3,main="Years of gainable life")
text(x=.82,y=0.74,"USPSTF")
text(x=.46,y=0.82,"Life gained based")
text(x=.40,y=0.17,"Risk-based")

text(x=.5,y=.51,"396k\n(33.4%)",cex=.75)
arrows(x0=.58, y0=.20, x1 = .61, y1 = .23, length = 0.025, angle = 30,
       code = 2)
text(x=.58,y=.18,"13k\n(1.1%)",cex=.75)
arrows(x0=.65, y0=.795, x1 = .6, y1 = .76, length = 0.025, angle = 30,
       code = 2)
text(x=.65,y=.8175,"31k\n(2.6%)",cex=.75)
text(x=.24,y=.545,"124k\n(10.5%)",cex=.75)
arrows(x0=.24, y0=.78, x1 = .32, y1 = .76, length = 0.025, angle = 30,
       code = 2)
text(x=.2,y=.78,"23k\n(2.0%)",cex=.75)
arrows(x0=.2, y0=.27, x1 = .24, y1 = .3, length = 0.025, angle = 30,
       code = 2)
text(x=.2,y=.25,"30k\n(2.6%)",cex=.75)
text(x=.81,y=.48,"64k\n(5.4%)",cex=.75)
text(x=.78,y=.18,"502k not gained\n(42.4%)",cex=.75)
dev.off()