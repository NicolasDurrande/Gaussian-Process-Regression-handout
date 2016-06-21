#library(DiceKriging)
library(tikzDevice)
library(MASS)
library(scatterplot3d)
library(lattice)
source("a_functions.R")

###############################################
### FIG gauss vec

K <- matrix(c(1,2,2,7),2)
Y <- mvrnorm(700,c(0,2),K)
tikz('ch1_gaussvec1.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
plot(Y,xlab='$Y_1$',ylab='$Y_2$',asp=1,col=rgb(0,0,0,.5),cex.axis=2,cex.lab=2)
dev.off()
tools::texi2dvi('ch1_gaussvec1.tex',pdf=T)

K <- matrix(c(1,0,0,1),2)
Y <- mvrnorm(1000,c(0,0),K)
tikz('ch1_gaussvec2.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
plot(Y,xlab='$Y_1$',ylab='$Y_2$',asp=1,col=rgb(0,0,0,.5),cex.axis=2,cex.lab=2)
dev.off()
tools::texi2dvi('ch1_gaussvec2.tex',pdf=T)

K <- matrix(c(1,0.4,0.8,0.4,0.5,0.3,0.8,0.3,1),3)
Y <- mvrnorm(1000,c(0,0,0),K)
tikz('ch1_gaussvec3.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(1.5,5.1,0,1.5))
scatterplot3d(Y,xlab='$Y_1$',ylab='$Y_2$',zlab='$Y_3$',grid=FALSE,asp=1,color=rgb(0,0,0,.5),angle=70,cex.axis=1.5,cex.lab=2)
dev.off()
tools::texi2dvi('ch1_gaussvec3.tex',pdf=T)

K <- matrix(c(4,-2,-2,1.5),2)
Y <- mvrnorm(1500,c(0,0),K)
for(i in 1:1500){
  if(runif(1)>.7 ) Y[i,1] <- -Y[i,1] 
}
tikz('ch1_gaussvec4.tex', standAlone = TRUE, width=5,height=5)
par(mar=c(4.5,5.1,1.5,1.5))
plot(Y,xlab='$Y_1$',ylab='$Y_2$',asp=1,col=rgb(0,0,0,.5),cex.axis=2,cex.lab=2)
dev.off()
tools::texi2dvi('ch1_gaussvec4.tex',pdf=T)

###############################################
### multivariate pdf

xgrid <- seq(-8.2,8.2,length.out=31)
Xgrid <- as.matrix(expand.grid(xgrid,xgrid))

K <- matrix(c(2,2,2,7),2)
mu <- matrix(c(0,0))
L <- GaussDensity(Xgrid,mu,K)

tikz('ch1_pdf1.tex', standAlone = TRUE, width=5,height=5)
par(mar=c(4.5,5.1,1.5,1.5))
persp(xgrid,xgrid,matrix(100*L,length(xgrid)),nticks =2,theta = 25, phi = 25,border=darkblue,ticktype="detailed",xlab='$x_1$',ylab='$x_2$',zlab='$100 \\times f_Y(x)$',cex.axis=1.5,cex.lab=2)
dev.off()
tools::texi2dvi('ch1_pdf1.tex',pdf=T)


xgrid <- seq(-5,5,length.out=101)
Xgrid <- as.matrix(expand.grid(xgrid,xgrid))
K <- matrix(c(2,2,2,7),2)
mu <- matrix(c(0,0))
L <- GaussDensity(Xgrid,mu,K)
sig = (2*pi)^(1/2)*det(K)
levels = dnorm(qnorm(seq(11/20,19/20,length.out=9),0,sig),0,sig)

tikz('ch1_pdf2.tex', standAlone = TRUE, width=5,height=5)
par(mar=c(4.5,5.1,1.5,1.5))
contour(xgrid,xgrid,matrix(L,length(xgrid)),drawlabels =FALSE,levels =levels,xlab='$x_1$',ylab='$x_2$',lwd=2,cex.axis=1.5,cex.lab=2,asp=1)
dev.off()
tools::texi2dvi('ch1_pdf2.tex',pdf=T)


###############################################
### 2D multivariate conditioning

xgrid <- seq(-8.2,8.2,length.out=31)
ygrid <- xgrid[10:31]
Xgrid <- as.matrix(expand.grid(xgrid,ygrid))

K <- matrix(c(7,5,5,7),2)
mu <- matrix(c(0,0))
L <- GaussDensity(Xgrid,mu,K)

x1 <- seq(-8.2,8.2,length.out=51)
x2 <- 0*x1 + ygrid[1]
Lc <- GaussDensity(cbind(x1,x2),mu,K)

muc <- K[1,2]/K[2,2]*ygrid[1]
sigmac <- sqrt(K[1,1] - K[1,2]/K[2,2]*K[2,1])
tikz('ch1_condpdf1.tex', standAlone = TRUE, width=5,height=5)
par(mar=c(.5,.5,.5,.5))
pmat <- persp(xgrid,ygrid,matrix(100*L,length(xgrid)),nticks =2,theta = 25, phi = 30,expand=.9, border=darkblue,xlab='$x_1$',ylab='$x_2$',zlab='$f_Y$',cex.axis=1.5,cex.lab=2,xlim=range(Xgrid),ylim=range(Xgrid),zlim=c(0,.7))
polygon(trans3d(x1, x2, 100*Lc, pmat),col=darkblue,border="black")
polygon(trans3d(8.2*c(-1,-1,1,1), rep(ygrid[1],4), c(0,.7,.7,0), pmat),border=NA,col=lightblue)
lines(trans3d(rep(muc,2),rep(ygrid[1],2),c(0,max(Lc)*100),pmat),lty=2,lwd=1.5)
lines(trans3d(c(muc,muc-sigmac),rep(ygrid[1],2),rep(0.18,2),pmat),lty=2,lwd=1.5)
text(trans3d(muc,ygrid[1]-1,0,pmat),'$\\mu_c$',cex=1.5)
text(trans3d(muc-1.5*sigmac,ygrid[1]-1,0.23,pmat),'$\\sqrt{\\Sigma_c}$',cex=1.5)
dev.off()
tools::texi2dvi('ch1_condpdf1.tex',pdf=T)

###############################################
### 3D multivariate conditioning

##
K <- matrix(c(1,0.4,0.8,0.4,1,0.3,0.8,0.3,1),3)
Y <- mvrnorm(3000,c(0,0,0),K)
Y <- Y[abs(Y[,1])<3,]
Y <- Y[abs(Y[,2])<3,]
Y <- Y[abs(Y[,3])<3,]
Y1 <- Y[Y[,3]<1.5,]
Y2 <- Y[Y[,3]>1.5,]

xgrid <- seq(-3,3,length.out=51)
Xgrid <- cbind(as.matrix(expand.grid(xgrid,xgrid)),2+0*xgrid)

mu <- matrix(c(0,0,0))
L <- GaussDensity(Xgrid,mu,K)

Kcond <- K[2:3,2:3] - matrix(K[2:3,1]) %*% K[1,2:3] / K[1,1]
sig = (2*pi)^(1/2)*det(Kcond)
levels = dnorm(qnorm(seq(10/20,19/20,length.out=9),0,sig),0,sig)

tikz('ch1_condpdf2.tex', standAlone = TRUE, width=5,height=5)
par(mar=c(.5,.5,.5,.5))
pmat <- persp(0:1, 0:1, matrix(0,2,2), xlim=c(-3,3), ylim=c(-3,3), zlim=c(-3,3),expand=.9, 
              theta=25, phi=30, xlab="$x_1$", ylab="$x_2$", zlab="$x_3$",border=NA,,cex.axis=1.5,cex.lab=2)
points(trans3d(Y1[,1],Y1[,2],Y1[,3], pmat),col= rgb(0,0,0,.5),pch=16, cex=.65)
polygon(trans3d(3*c(-1,-1,1,1), 3*c(-1,1,1,-1), rep(2,4), pmat),border=NA,col=lightblue)
cl <- contourLines(xgrid,xgrid,matrix(L,51)/max(L)*levels[1],levels=levels)
for(i in 1:length(cl)){
  polygon(trans3d(cl[[i]]$x, cl[[i]]$y, 2+0*cl[[i]]$x, pmat),border='black',col=rgb(32/255,74/255,135/255,.3))
} 
points(trans3d(Y2[,1],Y2[,2],Y2[,3], pmat),col= rgb(0,0,0,.5),pch=16, cex=.65)
dev.off()
tools::texi2dvi('ch1_condpdf2.tex',pdf=T)


###############################################
### FIG brown path
n <- 300
x <- matrix(seq(from=0, to=1, length=n))
K <- kBrown(x,x)
mp <- 0*x
vp <- diag(K)
y <- mvrnorm(5,mp,K)

tikz('ch1_GPpath0.tex', standAlone = TRUE, width=5, height=5)
plotGPR(x,mp,vp)
for (i in 1:5) {
  lines(x, y[i,], col=darkblue)
}
dev.off()
tools::texi2dvi('ch1_GPpath0.tex',pdf=T)

### FIG 2b Gaussian samples
n <- 100
x <- matrix(seq(from=0, to=1, length=n))

K <- kGauss(x,x)
mp <- 0*x
vp <- diag(K)
y <- mvrnorm(5,mp,K)

tikz('ch1_GPpath2.tex', standAlone = TRUE, width=5, height=5)
plotGPR(x,mp,vp)
for (i in 1:5) {
	lines(x, y[i,], col=darkblue)
}
dev.off()
tools::texi2dvi('ch1_GPpath2.tex',pdf=T)

### FIG 2c Trajectoires Matern 3/2
n <- 150
x <- matrix(seq(from=0, to=1, length=n))
K <- kMat32(x,x)
mp <- -2*x+3*x^2
vp <- diag(K)
y <- mvrnorm(5,mp,K)


tikz('ch1_GPpath3.tex', standAlone = TRUE, width=5, height=5)
plotGPR(x,mp,vp)
for (i in 1:5) {
  lines(x, y[i,], col=darkblue)
}
dev.off()
tools::texi2dvi('ch1_GPpath3.tex',pdf=T)

### FIG 2d 2D matern 5/2 samples
n <- 30
xgrid <- seq(from=0, to=1, length=n)
Xgrid <- as.matrix(expand.grid(xgrid,xgrid))

K <- kMat52(Xgrid,Xgrid,c(1,.5,.5))
mp <- matrix(0,n^2,1)
y <- mvrnorm(1,m,K)

tikz('ch1_GPpath4.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
persp(xgrid,xgrid,matrix(y,length(xgrid)),col=lightblue,border=darkblue,nticks =2,theta = 25, phi = 25,ticktype="detailed",xlab='$x_1$',ylab='$x_2$',zlab='$Z(x)$',cex.axis=1.5,cex.lab=2)
dev.off()
tools::texi2dvi('ch1_GPpath4.tex',pdf=T)
