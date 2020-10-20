########################################################
## Draw figure 3 for manuscript                       ##
## James Liley, 6 Oct 2020                            ##
########################################################
##
## Throughout, we have
## xs and xa are univariate
## xl is constant at 0, and gL(.,xl)=xl
## g denotes g^A_e
## f denotes f_e



########################################################
## Packages, scripts, auxiliary functions, output     ##
##  directories                                       ##
########################################################

# Auxiliary functions
logit=function(x,sc=1) 1/ (1+ exp(-sc*x))
ilogit=function(x,sc=1) -log(1/x  - 1)/sc

# Output directory
save_dir="./"

# Comment/uncomment to save 
#pdf(paste0(save_dir,"convergence_plot.pdf"),width=5,height=7.5)
jpeg(paste0(save_dir,"convergence_plot.jpg"),width=5,height=7,units="in",res=400)


########################################################
## Define functions and parameters                    ##
########################################################

# Define reasonable values for f and g

# Rationale: additive risks
f=function(xs,xa) logit(xs+xa)

# Rationale: we intervene by lowering xa when rho>0.5, but 
#  since we have limited resources for intervention, we 
#  must also allow xa to increase when rho<0.5. We assume
#  that we can intervene more effectively when xa is 
#  high; if xa is low already we cannot make it much lower.
g=function(rho,xa) {
  # xa + 5*(logit(xa))*(rho-0.5)
  (xa + 0.5*(xa+sqrt(1+xa^2)))*(1-rho) + (xa - 0.5*(xa+sqrt(1+xa^2)))*rho
}

# Tolerance. When 
#    delta_e:= |rho_e(xs,xa)-rho_{e-1}(xs,xa)|<tol
#  we say rho has converged. If 
#    ||delta_e|-|delta_{e-1}||<tol  and |delta_e|> tol2
#  we say rho has diverged.
tol=0.01
tol2=0.05

# Number of iterations to consider in total
emax=100




########################################################
## Compute number of model iterations to convergence  ##
##  or divergence                                     ##
########################################################

# Plot limits
# We are not very interested in samples with low 
#  actionable risk, so we have a shorter plot range
xslim=c(-8,8)
xalim=c(-2,8) 

# Resolution
res=800

# Evaluate at these values of xs, xa
xs1=seq(xslim[1],xslim[2],length=res)
xa1=seq(xalim[1],xalim[2],length=res)

# Vectors of co-ordinates to evaluate
xsv=as.vector(outer(xs1,rep(1,length(xa1))))
xav=as.vector(outer(rep(1,length(xs1)),xa1))

# This matrix will be set to the number of 
#  iterations
icount=matrix(0,res,res)

# This matrix will be set to an indiciator
#  of whether rho(xs,xa) converges or 
#  diverges
iconv=matrix(-1,res,res)

# P(Y|X(0)=(xs,xa)) at epoch 0
p0=f(xsv,xav)

# This matrix will be used to record differences
#  in rho between epochs
delta=rep(0,length(xav))


for (e in 1:emax) {
  p1=f(xsv,g(p0,xav)) # p1 is now P(Y|X(0)=(xs,xa)) at epoch e
  
  delta_e=p1-p0 # difference in P(Y|...) between epoch i and epoch e-1
  
  wconv=which(abs(delta_e)< tol) # convergence at these coordinates
  wdiv=which(abs(abs(delta_e)-abs(delta)) < tol & (abs(delta_e) > tol2)) # diverged and settled
  wx=unique(c(wconv,wdiv))
  
  # Update matrices
  icount[wx]=1+icount[wx] 
  iconv[wconv]=1
  iconv[wdiv]=-1
  
  # Reset
  p0=p1 
  delta=delta_e  
}
# Change icount matrix into number of iterations until settling
icount=emax-icount




########################################################
## Draw plot                                          ##
########################################################

## Set up panels
xdim=1
layout(rbind(c(1,2),c(3,3),c(3,3)),widths=c(xdim,xdim,2*xdim),heights=c(xdim,xdim,xdim))
par(mar=c(4,4,1,1))

## Draw a plot showing f
plot(0,type="n",xlim=c(-4,4),ylim=c(-4,4),xaxs="i",yaxs="i",
  xlab=expression(paste("Non-actionable (static) risk factor, X"^s)),
  ylab=expression(paste("Actionable risk factor, X"^a)))
xas=seq(-4,4,length=res)
xcol=colorRampPalette(c("red","blue"))(100)
image(xas,xas,outer(xas,xas,f),col=xcol,add=T)
legend("topright",c("0","0.5","1"),pch=16,col=c(xcol[1],xcol[50],xcol[100]),
  bg="white",box.col="white",
  title=expression(paste("f(x"^s,",x"^a,")")))

## Draw a plot showing g
plot(0,type="n",xlim=c(-4,4),ylim=c(-4,4),xaxs="i",yaxs="i",
  xlab=expression(paste("Pre-intervention X"^a,"(0)")),
  ylab=expression(paste("Post-intervention X"^s,"(1)")))
xas=seq(-4,4,length=res)
rvals=seq(0,1,length=5)
rcols=c("red","pink","black","lightblue","blue")
for (i in 1:length(xas)) lines(xas,g(rvals[i],xas),col=rcols[i])
legend("topleft",legend=rvals,pch=16,col=rcols,bty="n",
  title=expression(paste(rho," value")))

# Example, if needed
#lines(c(2,2),c(-4,g(0.75,2)),lty=2)
#lines(c(-4,2),c(g(0.75,2),g(0.75,2)),lty=2)
#text(2,-3.9,expression(paste("(X"[a],"(0),ρ)=(2,0.75) ")), adj=c(1,0))
#text(-4,g(0.75,2)+0.1,expression(paste(" X"[a],"(1)")),adj=c(0,0))




# Initialise main plot
plot(0,type="n",xlim=xslim,ylim=xalim,xaxs="i",yaxs="i",
  xlab=expression(paste("Non-actionable (static) risk factor, X"^s)),
  ylab=expression(paste("Actionable risk factor, X"^a)))

# Plot convergent and divergent regions separately
ilim_pos=10; ilim_neg=10
mat_pos=icount; mat_pos[which(iconv<0)]=NA
mat_neg=icount; mat_neg[which(iconv>0)]=NA
mat_pos[which(mat_pos>ilim_pos)]=ilim_pos
mat_neg[which(mat_neg> ilim_neg)]= ilim_neg

cspread=1
icol_pos=colorRampPalette(c("white",rep("red",cspread),rep("darkred",cspread)))(50)
icol_neg=colorRampPalette(c("black",rep("blue",cspread),rep("lightblue",cspread)))(50)
image(xs1,xa1,mat_pos,col=icol_pos,add=T)
image(xs1,xa1,mat_neg,col=icol_neg,add=T)


# Add legend
legend("bottomleft",c("Div. fast","Div. slowly","Conv. slowly", "Conv. fast"),title="Behaviour",
  col=c("black","blue","red","white"),pch=16,bg="lightgray",bty="n")

# Add gradient; kludge this
ncol=length(icol_pos)+length(icol_neg)
xleg=xslim[1]+c(0.3,0.4)
yleg=xalim[1]+seq(0.3,1.4,length=ncol)
zleg=rbind(1:ncol,1:ncol)
image(xleg,yleg,zleg,col=c(icol_pos,icol_neg),add=T)


########################################################
## Draw sub-plots                                     ##
########################################################

##' @name h
##' Updating function for rho_e
##' @param xs value of xs
##' @param xa value of xa
##' @param r values of rho_{e-1}
h=function(xs,xa,r) f(xs,g(r,xa))

##' @name re
##' Function to compute rho_e. 
##' @param xs test point xs value
##' @param xa test point xa value
##' @param e number of iterations
re=function(xs,xa,e) {
  if (e==0) return(f(xs,xa)) else return(h(xs,xa,re(xs,xa,e-1)))
}

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
  
  # Evaluate rho_e at (xs,xa) for epochs 1-10. Inefficient, but at this scale OK.
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
  lines(xtrans(xv),ytrans(h(xs,xa,xv)))
  lines(xtrans(c(-0.05,1)),ytrans(c(-0.05,1)))
  
  # Draw cobweb  
  x=f(xs,xa); y=0
  for (i in 1:25) {
    y2=h(xs,xa,x)
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

# Test points
tp1=c(-2,4)
tp2=c(0.5,3)
tp3=c(1,1)

# Draw small boxes
plotp(tp1[1],tp1[2],xslim[2]-1.1*scx,xalim[2]-0.1*scy,scx,scy)
plotp(tp2[1],tp2[2],xslim[2]-1.1*scx,mean(xalim) + 0.75*scy,scx,scy)
plotp(tp3[1],tp3[2],xslim[2]-1.1*scx,xalim[1]+1.6*scy,scx,scy)

# Links
points(tp1[1],tp1[2],col="white",pch=16,cex=0.5)
points(tp2[1],tp2[2],col="white",pch=16,cex=0.5)
points(tp3[1],tp3[2],col="black",pch=16,cex=0.5)

# brief function to draw line with fixed length gap at end
gline=function(x,y,sub=0.9,...) {
  dx=sqrt((x[2]-x[1])^2 + (y[2]-y[1])^2)
  lines(mean(x) + ((dx-sub)/dx)*(x-mean(x)),mean(y) + ((dx-sub)/dx)*(y-mean(y)),...)
}
gline(c(tp1[1],xslim[2]-scx),c(tp1[2],xalim[2]-0.75*scy),col="darkgray",lty=2)
gline(c(tp2[1],xslim[2]-scx),c(tp2[2],mean(xalim)),col="darkgray",lty=2)
gline(c(tp3[1],xslim[2]-scx),c(tp3[2],xalim[1]+0.75*scy),col="darkgray",lty=2)


# Complete
dev.off()
