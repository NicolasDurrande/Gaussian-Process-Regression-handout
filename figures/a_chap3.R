library(tikzDevice)
library(MASS)
source("a_newfunctions.R")


ftest <- function(x) sin(2*pi*x) + x

##################################################################"
### Optim Leave One Out
n <- 101
x <- matrix(seq(from=0, to=1, length=n))
y <- ftest(x)

X <- matrix(seq(0.1,0.9,length.out=5))
Y <- ftest(X)

plot(X,Y)

LooRss <- function(theta){
  n <- nrow(X)
  predLoo <- rep(0,n)
  for(i in 1:n)
     predLoo[i] <- GPR(X[i,,drop=FALSE],X[-i,,drop=FALSE],Y[-i,,drop=FALSE],kMat52,param=c(1,theta))$mean
  return(sum((predLoo-Y)^2))
}

opt = optim(0.2,LooRss)
K_1 <- solve(kMat52(X,X,param=c(1,opt$par)))
sig2opt <- 1/nrow(X) * t(Y) %*% K_1 %*% diag(1/diag(K_1)) %*% K_1 %*% Y
 
c(sig2opt,opt$par)

pred <- GPR(x,X,Y,kGauss,param=c(sig2opt,opt$par))

tikz('ch3_CV.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
plotGPR(x, pred, xlab='$x$',ylab='$Z(x)|Z(X)=F$')
points(X, Y, col='black', pch=4,lwd=3)
dev.off()
tools::texi2dvi('ch3_CV.tex',pdf=T,clean=TRUE)
file.remove('ch3_CV.tex')


##################################################################"
### Optim Maximum likelihood
library(kergp)
colnames(X) <- "x"
colnames(Y) <- "y"


covMat52 <- covMan(kernel = kMat52,
                      hasGrad = FALSE,
                      acceptMatrix = TRUE,
                      label = "Matern 5/2",
                      d = 1,
                      inputs = "x",
                      parLower = c(sigma2 = 0.0, theta = 0.0),
                      parUpper = c( sigma2 = Inf, theta = .5),
                      parNames = c("sigma2","theta"),
                      par = c(theta = 0.2, sigma2 = 1))

mygp <- gp(formula = y ~ 1, data = data.frame(Y, X),
           inputs = colnames(X),
           cov = covMat52,
           compGrad=FALSE,
           varNoiseLower = 1e-8, varNoiseUpper = 1e-5)

mygp$covariance@par

pred <- GPR(x,X,Y,kGauss,param=mygp$covariance@par)

tikz('ch3_MLE.tex', standAlone = TRUE, width=7, height=5)
par(mar=c(5,5,.75,.75),cex.axis=1.5,cex.lab=2)
plotGPR(x, pred, xlab='$x$',ylab='$Z(x)|Z(X)=F$')
points(X, Y, col='black', pch=4,lwd=3)
dev.off()
tools::texi2dvi('ch3_MLE.tex',pdf=T,clean=TRUE)
file.remove('ch3_MLE.tex')


