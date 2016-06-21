library(DiceKriging)
library(tikzDevice)
library(MASS)
source("a_functions.R")

##################################################################"
### Modeles de krigeage
n <- 5
m <- 101
x <- matrix(seq(from=0, to=1, length=m))
X <- matrix(seq(from=0.1, to=0.9, length=n))
F <- matrix(c(0.5, 0, 1.5, 3, 2))

model <- GPR(x,X,F,kGauss)
                
m <- model[[1]]
K <- model[[2]]
y <- mvrnorm(60,m,K)

tikz('ch2_GPR.tex', standAlone = TRUE, width=5, height=5)
plotGPR(x,m,K,"$Z(x)|Z(X)=F$")
plot_lim = par("usr")
points(X, F, pch=4, cex=1,lwd=3)
dev.off()
tools::texi2dvi('ch2_GPR.tex',pdf=T)

tikz('ch2_obs.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
plot(X, F, pch=4, cex=1,lwd=3,ylim=plot_lim[3:4],xlab="$x$",ylab="$F$",cex.axis=2,cex.lab=2)
dev.off()
tools::texi2pdf('ch2_obs.tex')

tikz('ch2_GPcond.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
plot(x, m, type="n", xlab="$x$",ylab="$Z(x)|Z(X)=F$", ylim=plot_lim[3:4],cex.axis=2,cex.lab=2)
points(X, F, pch=4, cex=1,lwd=3)
for (i in 1:60) {
  lines(x, y[i,], col=darkbluetr)
}
dev.off()
tools::texi2dvi('ch2_GPcond.tex',pdf=T)

##################################################################"
##################################################################"
### GPR = linear combination of basis functions GAUSS
n <- 5
m <- 101
x <- matrix(seq(from=-0.5, to=1.5, length=m))
X <- matrix(seq(from=0.1, to=0.9, length=n))
F <- matrix(c(0.5, 0, 1.5, 3, 2))

model <- GPR(x,X,F,kGauss)

m <- model[[1]]
K <- model[[2]]

tikz('ch2_GPRbasisfuncGauss.tex', standAlone = TRUE, width=7, height=5)
plotGPR(x,m,K,"$Z(x)|Z(X)=F$")
plot_lim = par("usr")
points(X, F, pch=4, cex=1,lwd=3)
dev.off()
tools::texi2dvi('ch2_GPRbasisfuncGauss.tex',pdf=T)

tikz('ch2_basisfuncGauss.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
matplot(x,kGauss(x,X),type="l",lty=1,ylim=range(F),col=darkblue,ylab="$k(x,X)$",cex.axis=1.5,cex.lab=2,lwd=2)
points(X, F, pch=4, cex=1,lwd=3)
matlines(t(matrix(rep(X,2),ncol=2)),range(F),lty=2,col='black')
dev.off()
tools::texi2dvi('ch2_basisfuncGauss.tex',pdf=T)


##################################################################"
##################################################################"
### GPR = linear combination of basis functions Brown
n <- 5
m <- 101
x <- matrix(seq(from=0, to=1., length=m))
X <- matrix(seq(from=0.1, to=0.9, length=n))
F <- matrix(c(0.5, 0, 1.5, 3, 2))

model <- GPR(x,X,F,kBrown)

m <- model[[1]]
K <- model[[2]]

tikz('ch2_GPRbasisfuncBrown.tex', standAlone = TRUE, width=7, height=5)
plotGPR(x,m,K,"$Z(x)|Z(X)=F$")
plot_lim = par("usr")
points(X, F, pch=4, cex=1,lwd=3)
dev.off()
tools::texi2dvi('ch2_GPRbasisfuncBrown.tex',pdf=T)

tikz('ch2_basisfuncBrown.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(4.5,5.1,1.5,1.5))
matplot(x,kBrown(x,X),type="l",lty=1,ylim=range(F),col=darkblue,ylab="$k(x,X)$",cex.axis=1.5,cex.lab=2,lwd=2)
points(X, F, pch=4, cex=1,lwd=3)
matlines(t(matrix(rep(X,2),ncol=2)),range(F),lty=2,col='black')
dev.off()
tools::texi2dvi('ch2_basisfuncBrown.tex',pdf=T)



##################################################################"
##################################################################"
### Modeles de krigeage avec bruit
n <- 5
m <- 101
x <- matrix(seq(from=0, to=1, length=m))
X <- matrix(seq(from=0.1, to=0.9, length=n))
F <- matrix(c(0.5, 0, 1.5, 3, 2))

model <- GPR(x,X,F,kGauss,paramNoise=0.001)

m <- model[[1]]
K <- model[[2]]

tikz('ch2_GPRnoise0001.tex', standAlone = TRUE, width=5, height=5)
plotGPR(x,m,K,"$Z(x)|Z(X)+N(X)=F$")
plot_lim = par("usr")
#points(X, F, pch=4, cex=1,lwd=3)
dev.off()
tools::texi2dvi('ch2_GPRnoise0001.tex',pdf=T)


##################################################################"
##################################################################"
### Multi-output GPR
n <- 5
m <- 101
xx <- matrix(seq(from=0, to=1, length=m))
x <- cbind(rbind(xx,xx),c(rep(0,m),rep(.05,m)))

X <- matrix(c(seq(from=0.1, to=0.9, length=n),seq(from=0.1, to=0.4, length=n)))
X <- cbind(X,c(rep(0,n),rep(.05,n)))
XA <- X[1:n,]
XB <- X[(n+1):(2*n),]

K <- kGauss(X,X)
Y <- mvrnorm(1,rep(0,nrow(X)),K)

model <- GPR(x,X,Y,kGauss)

mA <- model[[1]][1:m]
KA <- diag(model[[2]][1:m,1:m])
low95A <- mA - 1.96*sqrt(pmax(0,KA))
upp95A <- mA + 1.96*sqrt(pmax(0,KA))

mB <- model[[1]][(m+1):(2*m)]
KB <- diag(model[[2]][(m+1):(2*m),(m+1):(2*m)])
low95B <- mB - 1.96*sqrt(pmax(0,KB))
upp95B <- mB + 1.96*sqrt(pmax(0,KB))

tikz('ch2_multiGPR.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(.5,.5,.5,.5))
pmat <- persp(xx,xx,matrix(0,m,m),nticks =2,theta = 25, phi = 30,expand=.9, border=NA,col=NA,xlab='$t$',ylab='label',zlab='$T$',cex.axis=1.5,cex.lab=2,zlim=c(-2.5,2.5))

# model A
lines(trans3d(xx, .3+0*xx, mA, pmat),col=darkblue,lwd=2)
polygon(trans3d(c(xx,rev(xx)),rep(.3,2*m), c(upp95A,rev(low95A)), pmat),border=darkblue,col=lightblue)
points(trans3d(X[1:n,1],rep(.3,n),Y[1:n], pmat),pch=16)

# model B
lines(trans3d(xx, .6+0*xx, mB, pmat),col=darkblue,lwd=2)
polygon(trans3d(c(xx,rev(xx)),rep(.6,2*m), c(upp95B,rev(low95B)), pmat),border=darkblue,col=lightblue)
points(trans3d(XB[,1],rep(.6,n),Y[(n+1):(2*n)], pmat),pch=16)

dev.off()
tools::texi2dvi('ch2_multiGPR.tex',pdf=T)


##################################################################"
##################################################################"
### GPR with derivatives

### Modeles de krigeage avec bruit
m <- 101
x <- matrix(seq(from=0, to=1, length=m))
X <- matrix(c(.1,.3,.6))
Xd <- matrix(c(.3,.8))
F <- matrix(sin(2*pi*X)+X)
Fd <- matrix(1+cos(Xd))

Kpp <- kGauss(X,X)
Kpd <- dkGauss(X,Xd)
Kdd <- ddkGauss(Xd,Xd)


K <- rbind(cbind(Kpp,Kpd),cbind(t(Kpd),Kdd))
Kx <- cbind(kGauss(x,X),dkGauss(x,Xd))

m <- Kx %*% solve(K) %*% rbind(F,Fd)
v <- diag(kGauss(x,x) -Kx %*% solve(K) %*% t(Kx) )

md <- cbind(kGauss(Xd,X),dkGauss(Xd,Xd)) %*% solve(K) %*% rbind(F,Fd)
tikz('ch2_derGPR.tex', standAlone = TRUE, width=5, height=5)
plotGPR(x,m,v,"$Z(x)|Z(X)=F,Z'(X')=F'$")
lines(Xd[1]+0.05*c(-1,1),md[1]+Fd[1]*0.05*c(1,-1),col='red',lwd=2.5)
lines(Xd[2]+0.05*c(-1,1),md[2]+Fd[2]*0.05*c(1,-1),col='red',lwd=2.5)
points(X, F, pch=4, cex=1,lwd=3)
dev.off()
tools::texi2dvi('ch2_derGPR.tex',pdf=T)



