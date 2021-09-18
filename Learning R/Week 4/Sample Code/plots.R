decimal.ticks <- function(fn.minexp = -4, fn.maxexp = +4){


  yy = c()
  for (ii in (fn.minexp:fn.maxexp)){
    ymin <- 10^ii;
    yy <- c(yy, seq(ymin, 10*ymin, by = ymin) )
  }
  decimal.ticks <- unique(yy)

}


plot.oep <- function ( lim.x2, lim.y1, filename, legend.text, leg.cex, wp.poi, wp.neg, wp.mix, wp.rl8  ) {


   lim.x1 = 0
   ax.size = 2
   lab.size = 1.6
   png(filename,width=1000,height=1000,type='cairo1')
   idx = c(2:length(wp.poi$OEP))
   plot(wp.poi$xx,wp.poi$OEP[idx],log=c('y'),type='l', lty=2,lwd=4,col='blue',
        axes = F, bg = 'transparent', cex = 1.5, xlab = '', ylab = '', ylim=c(lim.y1,1),xlim=c(lim.x1,lim.x2))
   lines(wp.neg$xx,wp.neg$OEP[idx],col='red',lty=2,lwd=4)
   lines(wp.mix$xx,wp.mix$OEP[idx],col='black',lty=1,lwd=4)
   lines(wp.rl8$xx, wp.rl8$OEP[idx], col='gray',lty=1,lwd=4)
   
   legend('topright',legend.text,
          lwd=rep(4,4,4,4),lty=c(2,2,1,1),col=c('blue','red','black','gray'),
          cex=leg.cex)

   xticks <- seq(0,lim.x2, by = 1e9)
   lab = xticks/1e+9
   axis(side=1, at = xticks, labels = lab, las = 1, tck = 0.02, cex.axis = ax.size)
     #mtext(side = 1, 'Hazard', cex = lab.size, padj = +2)

   xticks.2 <- seq(0,lim.x2, by = 1e+8)
   abline(v = xticks.2,lwd=.2)

   axis(side=1, at = xticks.2, labels = NA, tck = 0.01)
   ## X Axis (top)
   axis(side=3, at = xticks, labels = NA, tck = 0.02)
   axis(side=3, at = xticks.2, labels = NA, tck = 0.01)


   ## EP Y Axis
   yticks <- 10^seq(-10, 0)
   axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.02, cex.axis = ax.size)
   mtext(side = 2, 'EF', cex = lab.size, padj = -3)

   yticks.2 <- decimal.ticks()
   abline(h = yticks.2,lwd=.2)

   axis(side=2, at = yticks.2, labels = NA, las = 1, tck = -0.01, cex.axis = ax.size)

   ## RP Y Axis
   my.rps <- c(2,5,10,20,25,50,100,200,250,500,1000,5000)
   yticks.lbl <- c(1,my.rps)
   yticks<- 1./yticks.lbl
   axis(side=4, at = yticks, labels = yticks.lbl, las = 1, tck = -0.01, pos=  lim.x2 , cex.axis = ax.size)
   mtext(side = 4, 'RP', cex = lab.size, padj = 3)
   xx.txt <- lim.x2/3
   yy.txt <- 1/400
   my.lab = paste(country)
   text(xx.txt, yy.txt, my.lab, cex = 3, pos = 2, col = 'black', bg = 'white' )
   dev.off()

}

plot.aep <- function ( lim.x2, lim.y1, filename, legend.text, leg.cex, wp.poi, wp.neg, wp.mix, wp.rl8  ) {


   lim.x1 = 0
   ax.size = 2
   lab.size = 1.6
   png(filename,width=1000,height=1000,type='cairo1')
   idx = c(2:length(wp.poi$AEP))
   plot(wp.poi$xx,wp.poi$AEP[idx],log=c('y'),type='l', lty=2,lwd=4,col='blue',
        axes = F, bg = 'transparent', cex = 1.5, xlab = '', ylab = '', ylim=c(lim.y1,1),xlim=c(lim.x1,lim.x2))
   lines(wp.neg$xx,wp.neg$AEP[idx],col='red',lty=2,lwd=4)
   lines(wp.mix$xx,wp.mix$AEP[idx],col='black',lty=1,lwd=4)
   lines(wp.rl8$xx, wp.rl8$AEP[idx], col='gray',lty=1,lwd=4)

   legend('topright',legend.text,
          lwd=rep(4,4,4,4),lty=c(2,2,1,1),col=c('blue','red','black','gray'),
          cex=leg.cex)

   xticks <- seq(0,lim.x2, by = 1e9)
   lab = xticks/1e+9
   axis(side=1, at = xticks, labels = lab, las = 1, tck = 0.02, cex.axis = ax.size)
     #mtext(side = 1, 'Hazard', cex = lab.size, padj = +2)

   xticks.2 <- seq(0,lim.x2, by = 1e+8)
   abline(v = xticks.2,lwd=.2)

   axis(side=1, at = xticks.2, labels = NA, tck = 0.01)
   ## X Axis (top)
   axis(side=3, at = xticks, labels = NA, tck = 0.02)
   axis(side=3, at = xticks.2, labels = NA, tck = 0.01)


   ## EP Y Axis
   yticks <- 10^seq(-10, 0)
   axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.02, cex.axis = ax.size)
   mtext(side = 2, 'EF', cex = lab.size, padj = -3)

   yticks.2 <- decimal.ticks()
   abline(h = yticks.2,lwd=.2)

   axis(side=2, at = yticks.2, labels = NA, las = 1, tck = -0.01, cex.axis = ax.size)

   ## RP Y Axis
   my.rps <- c(2,5,10,20,25,50,100,200,250,500,1000,5000)
   yticks.lbl <- c(1,my.rps)
   yticks<- 1./yticks.lbl
   axis(side=4, at = yticks, labels = yticks.lbl, las = 1, tck = -0.01, pos=  lim.x2 , cex.axis = ax.size)
   mtext(side = 4, 'RP', cex = lab.size, padj = 3)
   xx.txt <- lim.x2/3
   yy.txt <- 1/400
   my.lab = paste(country)
   text(xx.txt, yy.txt, my.lab, cex = 3, pos = 2, col = 'black', bg = 'white' )
   dev.off()

}

plot.check.aep <- function ( lim.x2, lim.y1, filename, legend.text, leg.cex, wp.poi, wp.neg, wp.mix, sol  ) {


   lim.x1 = 0
   ax.size = 2
   lab.size = 1.6
   png(filename,width=1000,height=1000,type='cairo1')
   idx = c(2:length(wp.poi$AEP))
   plot(wp.poi$xx,wp.poi$AEP[idx],log=c('y'),type='l', lty=1,lwd=4,col='blue',
        axes = F, bg = 'transparent', cex = 1.5, xlab = '', ylab = '', ylim=c(lim.y1,1),xlim=c(lim.x1,lim.x2))
   lines(wp.neg$xx,wp.neg$AEP[idx],col='red',lty=1,lwd=4)
   lines(wp.mix$xx,wp.mix$AEP[idx],col='black',lty=1,lwd=4)
   lines(sol$xx, sol$AEPpoi, col='blue',lty=2,lwd=3)
   lines(sol$xx, sol$AEPneg, col='red',lty=2,lwd=3)
   lines(sol$xx, sol$AEPmix, col='black',lty=2,lwd=3)


   legend('topright',legend.text,
          lwd=rep(4,4,4,3,3,3),lty=c(1,1,1,2,2,2),col=c('blue','red','black','blue','red','black'),
          cex=leg.cex)

   xticks <- seq(0,lim.x2, by = 1e9)
   lab = xticks/1e+9
   axis(side=1, at = xticks, labels = lab, las = 1, tck = 0.02, cex.axis = ax.size)
     #mtext(side = 1, 'Hazard', cex = lab.size, padj = +2)

   xticks.2 <- seq(0,lim.x2, by = 1e+8)
   abline(v = xticks.2,lwd=.2)

   axis(side=1, at = xticks.2, labels = NA, tck = 0.01)
   ## X Axis (top)
   axis(side=3, at = xticks, labels = NA, tck = 0.02)
   axis(side=3, at = xticks.2, labels = NA, tck = 0.01)


   ## EP Y Axis
   yticks <- 10^seq(-10, 0)
   axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.02, cex.axis = ax.size)
   mtext(side = 2, 'EF', cex = lab.size, padj = -3)

   yticks.2 <- decimal.ticks()
   abline(h = yticks.2,lwd=.2)

   axis(side=2, at = yticks.2, labels = NA, las = 1, tck = -0.01, cex.axis = ax.size)

   ## RP Y Axis
   my.rps <- c(2,5,10,20,25,50,100,200,250,500,1000,5000)
   yticks.lbl <- c(1,my.rps)
   yticks<- 1./yticks.lbl
   axis(side=4, at = yticks, labels = yticks.lbl, las = 1, tck = -0.01, pos=  lim.x2 , cex.axis = ax.size)
   mtext(side = 4, 'RP', cex = lab.size, padj = 3)
   xx.txt <- lim.x2/3
   yy.txt <- 1/400
   my.lab = paste(country)
   text(xx.txt, yy.txt, my.lab, cex = 3, pos = 2, col = 'black', bg = 'white' )
   dev.off()

}

