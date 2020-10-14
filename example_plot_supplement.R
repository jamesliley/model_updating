########################################################
## Draw `chaos' figures for AISTATS supplement        ##
## Author redacted, 6 Oct 2020                        ##
########################################################
##
## Throughout, we have
## xs and xa are univariate
## xl is constant at 0, and gL(.,xl)=xl
## g denotes g^A_e
## f denotes f_e

## Resolution
res=300

## Number of epochs
N=25

## Wrapper function; we will use two sets of values for f,g
plotfg=function(f,g,xslim,xalim,tp1=c(-4.8,4),tp2=c(-5.2,2),tp3=c(-7,1)) {
  
## Points to evaluate at
xs=seq(xslim[1],xslim[2],length.out=res)
xa=seq(xalim[1],xalim[2],length.out=res)
  
# ptn(n,xs,xa) = P(Y(n)|Xs(0,n)=xs,Xa(0,n)=xa)
# = probability of outcome in epoch e, given Xs and Xa at time 0
# We can't use recursion as number of calls will grow exponentially:
#  since rho(e+1)=f(g(rho,inv.f(rho))). We use dynamic programming instead.
ptn=function(n,xs,xa) {
  p0=f(xs,xa)
  p1=f(xs,g(p0,xa))
  if (n==0) return(p0)
  if (n==1) return(p1)
  if (n>1) {
    pn=outer(1:length(xs),rep(0,n+1)) # pn[,i+1]=risk at epoch i/score fitted at epoch i
    xan=outer(1:length(xa),rep(0,n+1)) # xan[,i+1]=xa at epoch i
    pn[,1]=p0 # score from epoch 1 is p0
    pn[,2]=p1 # score from epoch 1 is p0
    xan[,1]=xa # at epoch 0, xa stays as is
    xan[,2]=g(p0,xa) # at epoch 1, we intervene on p0
    for (i in 2:n) {
      xan[,i+1]=g(pn[,i],xan[,i]) # Intervene on previous xa on basis of previous risk score
      pn[,i+1]=f(xs,xan[,i+1]) # Fit new risk score
    }
    return(pn[,n+1])
  }
}


## Update function from rho_n to rho_(n+1). This requires inverting f.
# Inverse of f. Can only be done for one value of xs at a time
f_inv=function(x,xs=xs0) {
  xas=seq(-10,10,length=10000)
  y=f(xs,xas)
  return(approx(y,xas,x)$y)
}
h=function(x,xs=xs0,xa=xa0) f(xs,g(x,f_inv(x,xs=xs)))


##' @name re
##' Function to compute rho_e. 
##' @param xs test point xs value
##' @param xa test point xa value
##' @param e number of iterations
re=function(xs,xa,e) {
  if (e==0) return(f(xs,xa)) else return(h(re(xs,xa,e-1),xs,xa))
}

# Recall xs and xa are linear sequences across the plot range.
x0=outer(xs,xa,f)
xn=outer(xs,xa,function(x,y) ptn(N,x,y)) # matrix of values of true risk P(Y_N|X_N(0)=(xs,xa)) at epoch N. Equivalentl to predicted risk at epoch N+1.
xd=xn-0.5


plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Set risk",ylab="Actionable risk",
  xaxs="i",yaxs="i")
image(xs,xa,xd,col=icol2,breaks=seq(-max(abs(xd)),max(abs(xd)),length=1+length(icol2)),add=T)
#contour(xs,xa,xd,add=T)
legend("bottomleft",c("Positive","Negative"),
  col=c("red","blue"),pch=16,bg="lightgray")


##' @name plotp 
##' This function draws a cobweb plot and a zig-zag plot
##'  showing convergence or divergence of rho_e at some 
##'  point. 
##' @param xs test point xs value
##' @param xa test point xa value
##' @param locx top-right corner of plot
##' @param locy top-left corner of plot
##' @param scx x-scale for plot, on scale of original
##' @param scy y-scale for plot, on scale of original
##' @param ... passed to lines()
##' @returns null, draws plot (on top of existing)
plotp=function(xs,xa,locx,locy,scx=1,scy=1,...) {
  
  # Evaluate rho_e at (xs,xa) for epochs 1-20. Inefficient, but at this scale OK.
  rho_e=rep(0,20); 
  for (e in 1:length(rho_e)) rho_e[e]=re(xs,xa,e)
  
  # Some parameters
  ztop=1.1 # top margin of zigzag plot 
  zbottom=1.5
  
  # Add plot manually. Draw white square and rectangle
  polygon(c(locx,locx+scx,locx+scx,locx),c(locy,locy,locy-scy,locy-scy),
    col="lightgray",border=NA)
  polygon(c(locx,locx+scx,locx+scx,locx),
    c(locy-scy*ztop,locy-scy*ztop,locy-scy*zbottom,locy-scy*zbottom),
    col="lightgray",border=NA)
  
  # Transforms for x and y
  xtrans=function(x) locx + 0.05*scx + 0.95*scx*x
  ytrans=function(y) locy - scy + 0.05*scy + 0.95*scy*y
  
  # Draw x and y axes
  lines(xtrans(c(-0.05,1)),ytrans(c(0,0)))
  lines(xtrans(c(0,0)),ytrans(c(-0.05,1)))
  
  # Draw update function and xy line
  xv=seq(0,1,length=50)
  lines(xtrans(xv),ytrans(h(xv,xs,xa)))
  lines(xtrans(c(-0.05,1)),ytrans(c(-0.05,1)))
  
  # Draw cobweb  
  x=f(xs,xa); y=0
  for (i in 1:25) {
    y2=h(x,xs,xa)
    lines(xtrans(c(x,x)),ytrans(c(y,y2)),col="red")
    lines(xtrans(c(x,y2)),ytrans(c(y2,y2)),col="red")
    x=y2; y=x
  }
  
  # Draw plot underneath
  lines(locx + scx*seq(0,1,length=length(rho_e)),
    locy - scy*zbottom + 0.9*scy*(zbottom-ztop)*(rho_e-min(rho_e))/(max(rho_e)-min(rho_e)),
    col="red")
  
}

vscx=0.18 # size of small boxes
scx=vscx*(xslim[2]-xslim[1])
scy=vscx*(xalim[2]-xalim[1])


# Draw small boxes
plotp(tp1[1],tp1[2],xslim[2]-1.1*scx,xalim[2]-0.1*scy,scx,scy)
plotp(tp2[1],tp2[2],xslim[2]-1.1*scx,mean(xalim) + 0.75*scy,scx,scy)
plotp(tp3[1],tp3[2],xslim[2]-1.1*scx,xalim[1]+1.6*scy,scx,scy)

# Links
points(tp1[1],tp1[2],col="black",pch=16,cex=1)
points(tp2[1],tp2[2],col="black",pch=16,cex=1)
points(tp3[1],tp3[2],col="black",pch=16,cex=1)

# brief function to draw line with fixed length gap at end
gline=function(x,y,sub=0.9,...) {
  dx=sqrt((x[2]-x[1])^2 + (y[2]-y[1])^2)
  lines(mean(x) + ((dx-sub)/dx)*(x-mean(x)),mean(y) + ((dx-sub)/dx)*(y-mean(y)),...)
}
gline(c(tp1[1],xslim[2]-scx),c(tp1[2],xalim[2]-0.75*scy),col="black",lty=2)
gline(c(tp2[1],xslim[2]-scx),c(tp2[2],mean(xalim)),col="black",lty=2)
gline(c(tp3[1],xslim[2]-scx),c(tp3[2],xalim[1]+0.75*scy),col="black",lty=2)

}



### Plot 1
f=function(xs,xa,sc=1) logit(sc*(xs + xa))
g=function(r,xa) (xa + 0.5*(xa+sqrt(1+xa^2)))*(1-r) + (xa - 0.5*(xa+sqrt(1+xa^2)))*r


# Comment/uncomment to save 
save_dir="./"
pdf(paste0(save_dir,"chaos_plot_1.pdf"),width=5,height=5)
#jpeg(paste0(save_dir,"chaos_plot_1.jpg"),width=5,height=5,units="in",res=400)
plotfg(f,g,c(-8,2),c(-6,4),tp1=c(-1,3.5),tp2=c(-5.2,2),tp3=c(-6,0))
dev.off()



### Plot 2
f=function(xs,xa,sc=1) logit(sc*(xs + xa))
g=function(r,xa) (xa + 2*logit(xa))*(1-r) + (xa - 2*logit(xa))*r


# Comment/uncomment to save 
save_dir="./"
pdf(paste0(save_dir,"chaos_plot_2.pdf"),width=5,height=5)
#jpeg(paste0(save_dir,"chaos_plot_2.jpg"),width=5,height=5,units="in",res=400)
plotfg(f,g,c(-8,8),c(-8,8),tp1=c(-4.8,3.5),tp2=c(4,-1),tp3=c(3,-6))
dev.off()



