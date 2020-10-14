########################################################
## Code for web app to illustrate AISTATS manuscript  ##
## James Liley, 6 Oct 2020                            ##
########################################################
##
## The authors intend that should this manuscript be 
##  accepted, this web app will be published online as
##  a supplement to the manuscript.
##
## Throughout, we have
## xs and xa are univariate
## xl is constant at 0, and gL(.,xl)=xl
## g denotes g^A_e
## f denotes f_e


#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

## Some functions
logit=function(x,sc=1) 1/ (1+ exp(-sc*x))
ilogit=function(x,sc=1) -log(1/x  - 1)/sc

# Additional details
description_text="
This application ignores 'latent' 
risk factors (XL in manuscript).

The function governing true risk 
represents the process determining
the outcome Y. We model Y as being
determined at a given instant, 
following an intervention. 

The input should be interpretable 
as a function in R. Parameter xs is
'set' risk, xa is 'actionable' risk.

The functioning governing the 
intervention effect represents an
action taken on actionable risk 
factors. We assume that, in each 
epoch, the following events occur 
in order: 1 - predictive scores are 
calculated based on covariates; 2 - 
any interventions in response to 
this risk score are made, possibly
modifying actionable covariates; 3
- Y occurs, or does not occur, 
with probability determined by 
covariates after intervention; 4 -
a new predictive model is fitted to
observations of pre-intervention 
covariates and outcomes Y. Please
see manuscript for further detail.

The input should be interpretable as
a function in R. Parameter rho is 
predictive score, parameter xa is 
actionable risk"


library(shiny)

shinyServer(function(input, output) {
  
  # Update input values on button click
  erx=eventReactive(input$update, 
    list(inf=input$f,
      ing=input$g,
      ine=input$e,
      inxs=input$xs,
      inxa=input$xa,
      mode=input$mode,
      xl=input$xl)
  )
  
  # Show additional details if panel is opened
  output$descriptionText <- renderText(description_text)
  
  # Main plots
  output$distPlot <- renderPlot({
    
    # Get reactive input on button press
    ix=erx()
    xl=ix$xl
    
    # Assign functions f and g
    tf=paste0("f=function(xs,xa) ",ix$inf)
    tg=paste0("g=function(rho,xa) ",ix$ing)
    eval(parse(text=tf))
    eval(parse(text=tg))
    # now f and g are set to whatever was in the input
    
    # Mode: naive updating or successive adjuvancy
    dynamic.mode=ix$mode
    

    # Other parameters
    xs=seq(-xl,xl,length=100)
    xa=xs
    N=ix$ine
    xa0=ix$inxa; xs0=ix$inxs; 
    
    
    # ptn(n,xs,xa) = P(Y(n)|Xs(0,n)=xs,Xa(0,n)=xa)
    # = probability of outcome in epoch e, given Xs and Xa at time 0
    # Although the map from (xs,xa) after intervention is the same for
    #  naive updating and successive adjuvancy, this function takes 
    #  different forms, as the interventions is different.
    if (dynamic.mode=="Naive updating") {
      # Note recursion
      ptn=function(n,xs,xa) {
        if (n==0) return(f(xs,xa)) else return(f(xs,g(ptn(n-1,xs,xa),xa)))
      }
      
      ## Update function from rho_n to rho_(n+1)
      h=function(x,xs=xs0,xa=xa0) f(xs,g(x,xa))
      
    } else {
      
      # can't use recursion as number of calls will grow exponentially:
      #  since rho(e+1)=f(g(rho,inv.f(rho))). Use dynamic programming instead.
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
      
    }

    
    ## Set up plot panels    
    par(mfrow=c(5,2))

    # Recall xs and xa are linear sequences across the plot range.
    x0=outer(xs,xa,f)
    xn=outer(xs,xa,function(x,y) ptn(N,x,y)) # matrix of values of true risk P(Y_N|X_N(0)=(xs,xa)) at epoch N. Equivalentl to predicted risk at epoch N+1.
    xp=outer(xs,xa,function(x,y) ptn(N-1,x,y)) # matrix of values of true risk P(Y_N|X_N(0)=(xs,xa)) at epoch N-1. Equivalentl to predicted risk at epoch N.
    xz=outer(xs,xa,function(x,y) ptn(N+1,x,y)) # matrix of values of true risk P(Y_N|X_N(0)=(xs,xa)) at epoch N+1. Equivalentl to predicted risk at epoch N+2.
    xd=xn-xp
    xd2=xz-xn
    
    # Colours for plots
    icol=colorRampPalette(c("white","black"))(100)
    icol2=colorRampPalette(c("blue","white","red"))(100)
    
    
    ############################################################
    ## Plot 1. Illustrate function f(xs,xa)                   ##
    ############################################################
    ##
    ## Draw contour plot of f. Draw test point in red.
    ##
    plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Set risk",ylab="Actionable risk",
      main=paste0("P(Y) given covariates"))
    contour(xs,xa,outer(xs,xa,f),add=T)
    points(xs0,xa0,col="red",pch=16)

    
    
    ############################################################
    ## Plot 2. Illustrate function g(rho,xa)                  ##
    ############################################################
    ##
    ## Draw function g(rho,xa) for rho=0, rho=0.5, rho=1.
    ##  rho=0: change in xa if risk score is 0 (lowest risk)
    ##  rho=0.5; change in xa if risk score is 0.5 (equivocal risk)
    ##  rho=1; change in xa if risk score is 1 (highest risk)
    ##
    plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Pre-intervention actionable risk",
      ylab="Post-intervention actionable risk",
      main=paste0("Intervention effect"))
    lines(xa,g(0.5,xa),lty=1)
    lines(xa,g(0,xa),col="red")
    lines(xa,g(1,xa),col="blue")
    legend("bottomright",c(expression(paste(rho," = 0")),
      expression(paste(rho," = 0.5")),expression(paste(rho," = 1"))),
      col=c("red","black","blue"),lty=2)
    

    
    
    ############################################################
    ## Plot 3. Illustrate true risk at epoch N.               ##
    ############################################################
    ##
    ## Density plot of true risk at epoch N, or predicted risk 
    ##  at epoch N+1. Draw test point in red.
    ##
    plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Set risk",ylab="Actionable risk",
      main=paste0("True risk, epoch ",N))
    image(xs,xa,xn,col=icol,breaks=seq(0,1,length=1+length(icol)),add=T)
    points(xs0,xa0,col="red",pch=16)
    legend("bottomleft",c("0","1"),
      col=c("white","black"),pch=16,bg="lightgray")

    
    
    ############################################################
    ## Plot 4. Illustrate predicted risk at epoch N.          ##
    ############################################################
    ##
    ## Density plot of predicted risk at epoch N, or true risk 
    ##  at epoch N-1. Draw test point in red.
    ##
    plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Set risk",ylab="Actionable risk",
      main=paste0("Predicted risk, epoch ",N))
    image(xs,xa,xp,col=icol,breaks=seq(0,1,length=1+length(icol)),add=T)
    points(xs0,xa0,col="red",pch=16)
    legend("bottomleft",c("0","1"),
      col=c("white","black"),pch=16,bg="lightgray")
    

    
    ############################################################
    ## Plot 5. Illustrate predicted risk at epoch N.          ##
    ############################################################
    ##
    ## Density plot of predicted risk at epoch N, or true risk 
    ##  at epoch N-1. Draw test point in red.
    ##
    plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Set risk",ylab="Actionable risk",
      main=paste0("Diff. between predicted risk in epoch ",N-1," and epoch ",N))
    image(xs,xa,xd,col=icol2,breaks=seq(-max(abs(xd)),max(abs(xd)),length=1+length(icol2)),add=T)
    contour(xs,xa,xd,add=T)
    points(xs0,xa0,col="red",pch=16)
    legend("bottomleft",c("Positive","Negative"),
      col=c("red","blue"),pch=16,bg="lightgray")
    
    
    
    ############################################################
    ## Plot 6. Difference in predicted risk between epochs    ##
    ############################################################
    ##
    ## Density plot of difference in predicted risk at epoch N, 
    ##  and predicted risk at epoch N+1, or equivalently true 
    ##  risk at epoch N-1 and true risk at epoch N. Also draws 
    ##  a contour plot. Draws test point in red.
    ##
    plot(0,type="n",xlim=range(xs),ylim=range(xa),xlab="Set risk",ylab="Actionable risk",
      main=paste0("Diff. between predicted risk in epoch ",N," and epoch ",N+1))
    image(xs,xa,xd2,col=icol2,breaks=seq(-max(abs(xd2)),max(abs(xd2)),length=1+length(icol2)),add=T)
    contour(xs,xa,xd2,add=T)
    points(xs0,xa0,col="red",pch=16)
    legend("bottomleft",c("Positive","Negative"),
      col=c("red","blue"),pch=16,bg="lightgray")


    
    ############################################################
    ## Plot 7. Track risk scores for test point.              ##
    ############################################################
    ##
    ## Simple line plot of values of rho_e for test point
    ##
    xx=rep(0,25); 
    for (i in 1:25) xx[i]=ptn(i,xs0,xa0); 
    plot(xx,type="l",xlab="Epoch",ylab="True risk",main="Convergence of risk with epoch")


    
    
    ############################################################
    ## Plot 8. Cobweb plot for test point.                    ##
    ############################################################
    ##
    ## Simple cobweb plot of values of rho_e for test point. 
    ##  Can be an easier way to see limits of recursion. Note 
    ##  that non-XY line represents updating function for rho. 
    ##  The plot of the updating function may be inaccurate for 
    ##  successive adjuvancy.
    ##
    plot(0,type="n",xlim=0:1,ylim=0:1,xlab=expression(rho[n]),ylab=expression(rho[n+1]),
      main=paste0("Cobweb plot for recursion"))
    abline(h=0,lty=1); abline(v=0,lty=1)
    x01=seq(0,1,length=100)
    abline(0,1,col="red",lty=2); lines(x01,h(x01),col="blue",lty=2)
    x0=f(xs0,xa0); y0=0
    x=x0; y=y0
    for (i in 1:25) {
      lines(c(x,x),c(y,h(x)))
      lines(c(x,h(x)),c(h(x),h(x)))
      x=h(x); y=x
    }



    ############################################################
    ## Plot 9. Convergent/divergent areas                     ##
    ############################################################
    ##
    ## Plots regions of x-y plane for which rho_e converges or
    ##  diverges. Draws uncertain regions in gray. Draws test
    ##  point in red
    ##
    xs1=seq(-xl,xl,length=400); xa1=xs1
    xp=outer(xs1,xa1,function(x,y) 
      abs(ptn(25,x,y)-ptn(24,x,y)))
    mx=max(xp,na.rm=T)
    image(xa1,xs1,xp,col=c("white","gray","black"),breaks=c(-1,0.2*mx,0.8*mx,2),#xaxt="n",yaxt="n",
      xlab="Set risk",ylab="Actionable risk",main="Convergent region")
    legend("bottomleft",c("Convergent","Oscillatory","Divergent"),
      col=c("white","gray","black"),bg="lightgray",pch=16)
    points(xs0,xa0,col="red",pch=16)

    
    ############################################################
    ## Plot 10. Time to converge/diverge                       ##
    ############################################################
    ##
    ## This plot shows roughly the relative number of epochs 
    ##  for the behaviour of rho_e to settle; that is, converge
    ##  or diverge to a stable oscillation between 0 and 1. 
    ##  Test point in red.
    ##
    x0=0*outer(xs1,xa1);

    xsv=as.vector(outer(xs1,rep(1,length(xa1))))
    xav=as.vector(outer(rep(1,length(xs1)),xa1))

    p0=f(xsv,xav)
    
    xav0=xav
    for (i in 1:20) {
      if (dynamic.mode=="Naive updating") {
        p1=f(xsv,g(p0,xav))
      } else {
        xav=g(p0,xav0)
        p1=f(xsv,xav)
        xav0=xav
      }

      df=p1-p0
      p0=p1
      
      w=which(abs(0.5-abs(df))>0.45) # This will be true if difference is close to 0 or 1
      x0[w]=1+x0[w] # add 1 where this is true.
    }
    rx=20-range(x0)
    icol=colorRampPalette(c("black","blue","white"))(50)
    image(xa1,xs1,x0,col=icol,#xaxt="n",yaxt="n",
      xlab="Set risk",ylab="Actionable risk",main="Number of epochs to convergence or divergence")
    legend("bottomleft",c("Fewer","","More"),
      col=c("white","blue","black"),pch=16,bg="lightgray")
    points(xs0,xa0,col="red",pch=16)
  })
  
})
