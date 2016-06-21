library(tikzDevice)
library(MASS)
library(kergp)
library(colorspace)

source("a_newfunctions.R")

##################################################################"
### Finite dimensional kernels

kFD <- function(x,y,param=NULL){
  if(is.null(param)) param <- c(1,1,1/(2*pi))
  param[1]*outer(cos(x[,1]/param[3]),cos(y[,1]/param[3])) + param[1]*outer(sin(x[,1]/param[3]),sin(y[,1]/param[3])) + param[2]*outer(cos(.5*x[,1]/param[3]),cos(.5*y[,1]/param[3])) + param[2]*outer(sin(.5*x[,1]/param[3]),sin(.5*y[,1]/param[3])) 
} 

x <- matrix(seq(from=0, to=1, length=101))

## Prior
mp <- 0*x
K <- kFD(x,x)
prior <- list(mean=mp, cov=K)
y <- t(mvrnorm(15,mp,K))

tikz('ch5_finitedim0.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
plotGPR(x,prior,xlab='$x$',ylim=c(-3,3))
matplot(x, y, col=darkBlue, type='l',lty=1,add=TRUE)
dev.off()
tools::texi2dvi('ch5_finitedim0.tex',pdf=T,clean=TRUE)
file.remove('ch5_finitedim0.tex')

## Posterior 2
X <- matrix(c(0.2,0.6))
Y <- matrix(c(-1,1))
pred <- GPR(x,X,Y,kFD)
y <- t(mvrnorm(15,pred$mean,pred$cov))

tikz('ch5_finitedim2.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
plotGPR(x,pred,xlab='$x$',ylim=c(-3,3))
matplot(x, y, col=darkBlue, type='l',ylim=range(y),lty=1,add=TRUE)
points(X, Y, col='black', pch=4,lwd=3)
dev.off()
tools::texi2dvi('ch5_finitedim2.tex',pdf=T,clean=TRUE)
file.remove('ch5_finitedim2.tex')

## Posterior 4
X <- matrix(c(0.2,0.6,0.4,0.8))
Y <- matrix(c(-1,1,0,0))
pred <- GPR(x,X,Y,kFD)

tikz('ch5_finitedim4.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
plotGPR(x,pred,xlab='$x$',ylim=c(-3,3))
points(X, Y, col='black', pch=4,lwd=3)
dev.off()
tools::texi2dvi('ch5_finitedim4.tex',pdf=T,clean=TRUE)
file.remove('ch5_finitedim4.tex')

##################################################################"
### Tensor sum of kernels
ftest <- function(x) sin(4*pi*x[,1,drop=FALSE]) + cos(4*pi*x[,2,drop=FALSE]) + 2*x[,2,drop=FALSE]

kerSum <- covMan(kernel = kSumGauss,
               hasGrad = FALSE,
               acceptMatrix = TRUE,
               label = "sumGauss",
               d = 2,
               parLower = c(sigma21 = 1e-5, theta1 = 1e-5, sigma22 = 1e-5, theta2 = 1e-5),
               parUpper = c(sigma21 = Inf, theta1 = 2, sigma22 = Inf, theta2 = 2),
               parNames = c("sigma21","theta1","sigma22","theta2"),
               par = c(sigma21 = 1, theta1 = .2, sigma22 = 1, theta2 = .2))


xgrid <- matrix(seq(from=0, to=1, length=50))
Xgrid <- expand.grid(xgrid,xgrid)
Ygrid <- ftest(Xgrid)

x <- matrix(runif(10000),ncol=2)
X <- data.frame(kmeans(x,12)$centers)
X <- 1.1*(X-.5) + .5
Y <- ftest(X)

names(Xgrid) <- paste("x",1:2,sep='')
names(X) <- paste("x",1:2,sep='')
colnames(Y) <- 'y'

## plot test function
tikz('ch5_addfunc.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
image(x = xgrid, y = xgrid, z = matrix(ftest(as.matrix(Xgrid)), nrow = length(xgrid), ncol = length(xgrid)),
      col = diverge_hcl(100), xlab='$x_1$',ylab='$x_2$',zlim=c(-4,4))
contour(x = xgrid, y = xgrid, z = matrix(ftest(as.matrix(Xgrid)), nrow = length(xgrid), ncol = length(xgrid)),
        nlevels = 10,add=TRUE,drawlabels=FALSE,labcex=2,lwd=2)
points(x2 ~ x1, data = X, type = "p", pch = 16, col = "black", cex = 1.5)
dev.off()
tools::texi2dvi('ch5_addfunc.tex',pdf=T,clean=TRUE)
file.remove('ch5_addfunc.tex')

## Build model
addgp <- gp(formula = y ~ 1, data = data.frame(Y, X),
            inputs = names(X),
            cov = kerSum,
            compGrad=FALSE,
            varNoiseLower = 1e-4, varNoiseUpper = 1e-1)

summary(addgp)

## predict and plot
predaddgp <- predict(object = addgp, newdata = Xgrid, type = "SK")

tikz('ch5_addmean.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
image(x = xgrid, y = xgrid, z = matrix(predaddgp$mean, nrow = length(xgrid), ncol = length(xgrid)),
      col = diverge_hcl(100), xlab='$x_1$',ylab='$x_2$',zlim=c(-4,4))
contour(x = xgrid, y = xgrid, z = matrix(predaddgp$mean, nrow = length(xgrid), ncol = length(xgrid)),
      nlevels = 10,add=TRUE,drawlabels=FALSE,labcex=2,lwd=2)
points(x2 ~ x1, data = X, type = "p", pch = 16, col = "black", cex = 1.5)
dev.off()
tools::texi2dvi('ch5_addmean.tex',pdf=T,clean=TRUE)
file.remove('ch5_addmean.tex')

tikz('ch5_addsd.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
image(x = xgrid, y = xgrid, z = matrix(4*predaddgp$sd, nrow = length(xgrid), ncol = length(xgrid)),
      col = diverge_hcl(100), xlab='$x_1$',ylab='$x_2$',zlim=c(-4,4))
contour(x = xgrid, y = xgrid, z = matrix(predaddgp$sd, nrow = length(xgrid), ncol = length(xgrid)),
        nlevels = 5,add=TRUE,drawlabels=FALSE,labcex=2,lwd=2)
points(x2 ~ x1, data = X, type = "p", pch = 16, col = "black", cex = 1.5)
dev.off()
tools::texi2dvi('ch5_addsd.tex',pdf=T,clean=TRUE)
file.remove('ch5_addsd.tex')

### Legend stripe
tikz('ch5_addlegend.tex', standAlone = TRUE, width=0.75, height=5)
par(mar=c(5,.75,.75,2.5),cex.axis=1.5,cex.lab=2)
image(x = 1, y = seq(-4,4,length.out=101), z = t(matrix(seq(-4,4,length.out=101))),
      col = diverge_hcl(100), xlab='',ylab='',xaxt='n',yaxt='n')
axis(4)
dev.off()
tools::texi2dvi('ch5_addlegend.tex',pdf=T,clean=TRUE)
file.remove('ch5_addlegend.tex')

##################################################################"
### Sub-models of additive kernels

par1 <- addgp$covariance@par[1:2]
par2 <- addgp$covariance@par[3:4]

subm1 <- GPR(cbind(xgrid,1),X,as.matrix(Y),kSumGauss,param=c(par1,0,1),kernNoise=kSumGauss,paramNoise=c(0,1,par2))
subm2 <- GPR(cbind(xgrid,xgrid),X,as.matrix(Y),kSumGauss,param=c(0,1,par2),kernNoise=kSumGauss,paramNoise=c(par1,0,1))

tikz('ch5_addSubModel1.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
plotGPR(xgrid, subm1, xlab='$x_1$',ylab='')
points(X[,1], Y$y, col='black', pch=4,lwd=3)
dev.off()
tools::texi2dvi('ch5_addSubModel1.tex',pdf=T,clean=TRUE)
file.remove('ch5_addSubModel1.tex')

tikz('ch5_addSubModel2.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
plotGPR(xgrid, subm2, xlab='$x_2$',ylab='')
points(X[,2], Y$y, col='black', pch=4,lwd=3)
dev.off()
tools::texi2dvi('ch5_addSubModel2.tex',pdf=T,clean=TRUE)
file.remove('ch5_addSubModel2.tex')


##################################################################"
### Product of kernels
x <- matrix(seq(-1,1,0.01))
K <- kCos(x,x,c(1,0.07)) * kMat52(x,x,c(1,.3))
y <- t(mvrnorm(6,0*x,K))

tikz('ch5_prodKern.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
plot(x,K[,101],type='l',xlab="$x$",ylab="",lwd=3,col=darkBlue)
dev.off()
tools::texi2dvi('ch5_prodKern.tex',pdf=T,clean=TRUE)
file.remove('ch5_prodKern.tex')

tikz('ch5_prodSimu.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
matplot(x, y, col=sample(diverge_hcl(ncol(y))), type='l',ylim=range(y),xlab='$x$',ylab='',lty=1,lwd=1.5)
dev.off()
tools::texi2dvi('ch5_prodSimu.tex',pdf=T,clean=TRUE)
file.remove('ch5_prodSimu.tex')


##################################################################"
### kernels Rescaling

x <- matrix(seq(.2,8,0.01))
K <- kMat32(x,x)
y <- t(mvrnorm(6,0*x,K))

tikz('ch5_rescKern1.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
matplot(1/x, y, col=sample(diverge_hcl(ncol(y))), type='l',ylim=range(y),xlab='$x$',ylab='',lty=1,lwd=1.5)
dev.off()
tools::texi2dvi('ch5_rescKern1.tex',pdf=T,clean=TRUE)
file.remove('ch5_rescKern1.tex')

tikz('ch5_rescKern2.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.5,.5),cex.axis=1.5,cex.lab=2)
matplot(x, matrix(1/x^2,nrow(y),ncol(y))*y, col=sample(diverge_hcl(ncol(y))), type='l',ylim=range(y),xlab='$x$',ylab='',lty=1,lwd=1.5)
dev.off()
tools::texi2dvi('ch5_rescKern2.tex',pdf=T,clean=TRUE)
file.remove('ch5_rescKern2.tex')

##################################################################"
### linear transformation kernels

cppFunction('
        NumericVector kernFun(NumericVector x1, NumericVector x2, 
                                 NumericVector par){
        int n1 = x1.size();
        double S, d1, d2; 
        NumericVector K(1), h(n1);
        h = (abs(x1) - abs(x2)) / par[0];  // sugar function "abs"
        S = sum(h * h);                    // sugar "*" and "sum" 
        d2 = exp(-S);
        K[0] = par[1] * d2;
        d1 = 2 * K[0] * S / par[0];   
        K.attr("gradient") = NumericVector::create(Named("theta", d1),
                                                   Named("sigma2", d2));
        return K;
     }')

covSymGauss <- covMan(kernel = kernFun,
                      hasGrad = TRUE,
                      label = "argumentwise invariant",
                      d = 2,
                      parLower = c(theta = 0.0, sigma2 = 0.0),
                      parUpper = c(theta = Inf, sigma2 = Inf),
                      parNames = c("theta", "sigma2"),
                      par = c(theta = 0.5, sigma2 = 2))

nGrid <- 50; n <- nGrid^2; d <- 2
xGrid <- seq(from = -1, to = 1, length.out = nGrid)
Xgrid <- expand.grid(x1 = xGrid, x2 = xGrid)

Kmat <- covMat(object = covSymGauss, X = Xgrid,
               compGrad = FALSE, index = 1L)

library(MASS)
set.seed(1)
ygrid <- mvrnorm(mu = rep(0, n), Sigma = Kmat)

## sample
tikz('ch5_symSample.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
image(x = xGrid, y = xGrid, z = matrix(ygrid, nrow = nGrid, ncol = nGrid),
      col = diverge_hcl(100), xlab='$x_1$',ylab='$x_2$',zlim=c(-max(abs(ygrid)),max(abs(ygrid))))
contour(x = xGrid, y = xGrid, z = matrix(ygrid, nrow = nGrid, ncol = nGrid),
        nlevels = 10,add=TRUE,drawlabels=FALSE,labcex=2,lwd=2)
#points(x2 ~ x1, data = X, type = "p", pch = 16, col = "black", cex = 1.5)
dev.off()
tools::texi2dvi('ch5_symSample.tex',pdf=T,clean=TRUE)
file.remove('ch5_symSample.tex')

## Fit the Gaussian process model
nDesign <- 20
rowIndex <- sample(1:n)[1:nDesign]
X <- Xgrid[rowIndex,]
KXX <- covMat(object = covSymGauss, X = X,compGrad = FALSE, index = 1L)

set.seed(2)
y <- mvrnorm(mu = rep(0, nDesign), Sigma = KXX)

symgp <- gp(formula = y ~ 1, data = data.frame(y, X),
            inputs = names(X),
            cov = covSymGauss,
            parCovIni = c(0.1, 2),
            varNoiseLower = 1e-8, varNoiseUpper = 1e-8)

summary(symgp)

## -- predict and compare --
predSymgp <- predict(object = symgp, newdata = Xgrid, type = "SK")

tikz('ch5_symmean.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
image(x = xGrid, y = xGrid, z = matrix(predSymgp$mean, nrow = nGrid, ncol = nGrid),
      col = diverge_hcl(100), xlab='$x_1$',ylab='$x_2$',zlim=c(-max(abs(predSymgp$mean)),max(abs(predSymgp$mean))))
contour(x = xGrid, y = xGrid, z = matrix(predSymgp$mean, nrow = nGrid, ncol = nGrid),
        nlevels = 10,add=TRUE,drawlabels=FALSE,labcex=2,lwd=2)
points(x2 ~ x1, data = X, type = "p", pch = 16, col = "black", cex = 1.5)
dev.off()
tools::texi2dvi('ch5_symmean.tex',pdf=T,clean=TRUE)
file.remove('ch5_symmean.tex')

tikz('ch5_addsd.tex', standAlone = TRUE, width=5, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
image(x = xGrid, y = xGrid, z = matrix(predSymgp$sd, nrow = nGrid, ncol = nGrid),
      col = diverge_hcl(100), xlab='$x_1$',ylab='$x_2$',zlim=c(-max(abs(predSymgp$sd)),max(abs(predSymgp$sd))))
contour(x = xGrid, y = xGrid, z = matrix(predSymgp$sd, nrow = nGrid, ncol = nGrid),
        nlevels = 10,add=TRUE,drawlabels=FALSE,labcex=2,lwd=2)
points(x2 ~ x1, data = X, type = "p", pch = 16, col = "black", cex = 1.5)
dev.off()
tools::texi2dvi('ch5_addsd.tex',pdf=T,clean=TRUE)
file.remove('ch5_addsd.tex')
