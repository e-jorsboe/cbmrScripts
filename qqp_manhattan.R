qqpemil<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",...){
    x<-x[!is.na(x)]
    maxLogP<-25
    if(!missing(maxLogP)){
        x[x<10^-maxLogP]<-10^-maxLogP
    }
    N<-length(x)
    chi1<-qchisq(1-x,1)
    x<-sort(x)
    lambda<-round(median(chi1)/qchisq(0.5,1),2)
    e<- -log((1:N-0.5)/N,10)
    if(add){
        points(e,-log(x,10),...)
    } else{
        par(mar=c(2.5,2.5,1,0.5), mgp=c(1.6,0.6,0), font.main=1, cex.main=1)
        plot(e,-log(x,10),ylab=ylab,xlab=xlab,bty="n",...)
        abline(h=25, lty=3)
        text(-0.3,25.5,'Capped',adj=c(0,0), font=3)
        box(bty='L')
        abline(0,1,col=2,lwd=2)
    }

    ##title(paste("lambda=",lambda), cex=1.5)
    mtext(paste0("lamdbda = ",lambda),cex=2)

    if(ci){
        c95<-qbeta(0.95,1:N,N-(1:N)+1)
        c05<-qbeta(0.05,1:N,N-(1:N)+1)
        lines(e,-log(c95,10))
        lines(e,-log(c05,10))
    }
}


qqpemil999quan<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",...){
    x<-x[!is.na(x)]
    maxLogP<-25
    if(!missing(maxLogP)){
        x[x<10^-maxLogP]<-10^-maxLogP
    }
    N<-length(x)
    chi1<-qchisq(1-x,1)
    x<-sort(x)
    lambda<-round(quantile(chi1,probs=0.999)/qchisq(0.999,1),2)
    e<- -log((1:N-0.5)/N,10)
    if(add){
        points(e,-log(x,10),...)
    } else{
        par(mar=c(2.5,2.5,1,0.5), mgp=c(1.6,0.6,0), font.main=1, cex.main=1)
        plot(e,-log(x,10),ylab=ylab,xlab=xlab,bty="n",...)
        abline(h=25, lty=3)
        text(-0.3,25.5,'Capped',adj=c(0,0), font=3)
        box(bty='L')
        abline(0,1,col=2,lwd=2)
    }

    ##title(paste("lambda=",lambda), cex=1.5)
    mtext(paste0("lamdbda (99.9 % quantile) = ",lambda),cex=2)

    if(ci){
        c95<-qbeta(0.95,1:N,N-(1:N)+1)
        c05<-qbeta(0.05,1:N,N-(1:N)+1)
        lines(e,-log(c95,10))
        lines(e,-log(c05,10))
    }
}



manPlot<-function(x,chr,pass,sign,main,collar=c("darkblue","#67a9cf","orange","grey")){
  keep<-!is.na(x)
  x<-x[keep]
  chr<-chr[keep]
  x<-ifelse(x<1e-30,1e-30,x)
  col<-chr%%2+1
  pch<-rep(16,length(x))
  if(!missing(pass)){ ## NA
    col[!pass[keep]]<-4
    pch[!pass[keep]]<-1
  }
  par(mar=c(5.1, 5.1, 4.1, 1.1))

  maxy = max((-log(x,base=10)))
  if(maxy<8){
    maxy<-8
  }
  plot(-log(x,base=10),col=collar[col],ylab=expression(-log[10](italic(P))),xlab="Chromosomes",main=main,cex=1,lwd=2,pch=pch,axes=F,cex.lab=2,cex.main=2,ylim=c(0,maxy+0.2*maxy))
  box()
  axis(2,las=1,cex.axis=1.8)
  t<--log(0.05/length(x),base=10)
  print(0.05/length(x))
  abline(h=sign,lty=2,col=1)
  mtext(unique(chr),side=1,at=as.vector(tapply(1:length(chr),chr,mean)))

  ##          legend(0,t-0.5,"Bonferroni correction",lty=2,bty="n",col=1,cex=2)
  ## legend("topright",c("known","novel"),pch=c(4,1),bty="n",cex=2)

}

qqPlot<-function(x,pass,phe,main="IHIT"){

  par(mar=c(2.5,2.5,1,0.5), mgp=c(1.6,0.6,0), font.main=1, cex.main=1)
  keep<-!is.na(x)
  x<-x[keep]
  x<-ifelse(x<1e-30,1e-30,x)
  pass=pass[keep]
  maxy = max((-log(x,base=10)))

  qqpemil(x,pch=16,col=ifelse(pass,"darkblue","grey"),main=phe,las=1,cex.lab=2,cex.main=2,ylim=c(0,maxy+0.2*maxy),
      xlab=expression(Expected~~-log[10](italic(P))),ylab=expression(Observed~~-log[10](italic(P))),cex.axis=1.5)

}

qqPlot999quan<-function(x,pass,phe,main="IHIT"){

  par(mar=c(2.5,2.5,1,0.5), mgp=c(1.6,0.6,0), font.main=1, cex.main=1)
  keep<-!is.na(x)
  x<-x[keep]
  x<-ifelse(x<1e-30,1e-30,x)
  pass=pass[keep]
  maxy = max((-log(x,base=10)))

  qqpemil999quan(x,pch=16,col=ifelse(pass,"darkblue","grey"),main=phe,las=1,cex.lab=2,cex.main=2,ylim=c(0,maxy+0.2*maxy),
      xlab=expression(Expected~~-log[10](italic(P))),ylab=expression(Observed~~-log[10](italic(P))),cex.axis=1.5)

}
