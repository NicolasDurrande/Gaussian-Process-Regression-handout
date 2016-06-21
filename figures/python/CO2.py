import pylab as pb
import numpy as np
#import sys
#sys.path.insert(4,'/home/nicolas/Documents/ARTICLES/SheffieldML/HarmonicRKHS/GPy-0.2')
import GPy

np.random.seed(1)
pb.ion()

pb.close('all')

#################################################
# load dataset

DATA = np.loadtxt('CO2-DATA.txt',skiprows=1)
nrow = DATA[:,4].shape[0]
YEAR = DATA[:,0]
MONTH = np.linspace(1,nrow,nrow)
CO2 = DATA[:,4]

pb.figure()
X = (1958 + MONTH/12.)[:,None]
Y = CO2[:,None]
pb.plot(X,Y,'g')
pb.xlim((min(X),max(X)+20))
pb.ylim((310,440))
pb.savefig('CO2-data.pdf',bbox_inches='tight')

Xnew = 2013+ np.asarray(range(12*20))[:,None]/12.

##############################################
## create GP model for periodic and aperiodic part

kb = GPy.kern.bias(input_dim=1,variance=10.)
kl = GPy.kern.linear(input_dim=1,variances=10.)
kq = GPy.kern.linear(input_dim=1,variances=1.)*GPy.kern.linear(input_dim=1,variances=1e-4)
kg = GPy.kern.rbf(input_dim=1,variance=80.,lengthscale=20.)
kp = GPy.kern.periodic_Matern32(input_dim=1,variance=50.,lengthscale=.5,period=1.,n_freq=15,lower=Xnew[0,0],upper=Xnew[-1,0])
ka = GPy.kern.rbf(input_dim=1,variance=1.,lengthscale=.5)
kn = GPy.kern.white(input_dim=1,variance=20.)

###########################################
###########################################
## CO2 rbf

##############
indo = range(0,100,2) + range(109,660,5)
Xo = X[indo,:]
Yo = Y[indo,:]

##############
m = GPy.models.GPRegression(Xo,Yo,kernel = kg)
m['.*len'] = 60.
m['.*var'] = 1000000.
m['.*noise'] = .1

pb.figure()
Yp,Vp,lp,up = m.predict(Xnew)
GPy.util.plot.gpplot(Xnew,Yp,lp,up)
pb.plot(X,Y,'g')
pb.savefig('CO2-rbfb.pdf',bbox_inches='tight')

##############
m = GPy.models.GPRegression(Xo,Yo,kernel = ka)
m['.*len'] = .5
m['.*var'] = 50000.
m['.*noise'] = 1.

pb.figure()
Yp,Vp,lp,up = m.predict(Xnew)
GPy.util.plot.gpplot(Xnew,Yp,lp,up)
pb.plot(X,Y,'g')
pb.savefig('CO2-rbfa.pdf',bbox_inches='tight')


##############
m = GPy.models.GPRegression(Xo,Yo,kernel = kg + ka)
m['rbf_2_variance'] = 10000
m['rbf_1_variance'] = 5.
#m['.*noise'] = 1.

pb.figure()
Yp,Vp,lp,up = m.predict(Xnew)
GPy.util.plot.gpplot(Xnew,Yp,lp,up)
pb.plot(X,Y,'g')
pb.savefig('CO2-rbfab.pdf',bbox_inches='tight')

pb.figure()
m.kern.plot(2020,plot_limits=(2013,2033))

##############
m = GPy.models.GPRegression(Xo,Yo,kernel = kg + ka + kp + kq)
#m['rbf_2_variance'] = 10000
m['rbf_1_variance'] = .5
#m['.*<times>'] = 1e-4

pb.figure()
Yp,Vp,lp,up = m.predict(Xnew)
GPy.util.plot.gpplot(Xnew,Yp,lp,up)
pb.plot(X,Y,'g')
pb.savefig('CO2-rbfabpq.pdf',bbox_inches='tight')



# foo


# ###########################################
# ###########################################
# ## CO2 top optim param

# k = kb + kl  + kg + kp + ka + kq
# #k.constrain_fixed([-1])
# #k.tie_params('.*<times>')

# ##############
# indo = range(0,100,2) + range(110,660,10)
# Xo = X[indo,:]
# Yo = Y[indo,:]

# mo = GPy.models.GPRegression(Xo,Yo,kernel = k)

# mo.ensure_default_constraints()
# mo.constrain_fixed('.*_period')
# mo.optimize(messages=1,max_f_eval=500)

# ##############
# indo = range(0,100,1) + range(105,660,5)
# Xo = X[indo,:]
# Yo = Y[indo,:]

# mo2 = GPy.models.GPRegression(Xo,Yo,kernel = k)
# mo2[''] = mo['']

# mo2.ensure_default_constraints()
# mo2.constrain_fixed('.*_period')
# mo2.optimize(messages=1,max_f_eval=100)

# ##############
# m = GPy.models.GPRegression(X,Y,kernel = k)
# m[''] = mo2['']

# Yp,Vp,lp,up = m.predict(Xnew)
# GPy.util.plot.gpplot(Xnew,Yp,lp,up)
# pb.plot(X,Y,'g')
# pb.savefig('CO2-top.pdf',bbox_inches='tight')

# foo

# ##############################################
# # plots

# pb.figure(figsize=(4,4))
# pb.plot(Xnew,Ynew+mean,'r-', label='line 1', linewidth=1.5)
# m, tmp, lower, upper = mm.predict(Xnew)
# GPy.util.plot.gpplot(Xnew,m+mean, lower+mean, upper+mean)
# pb.plot([X[-1],X[-1]],[-10+mean,15+mean],'k:')
# year_ticks()
# pb.savefig('CO2a-t.pdf')

# pb.figure(figsize=(4,4))
# m, tmp, lower, upper = mm.predict(Xnew, slices=(True,False,False,False))
# GPy.util.plot.gpplot(Xnew,m+mean, lower+mean, upper+mean)
# pb.plot([X[-1],X[-1]],[-10+mean,15+mean],'k:')
# year_ticks()
# pb.savefig('CO2a-mp.pdf')

# pb.figure(figsize=(4,4))
# m, tmp, lower, upper = mm.predict(Xnew, slices=(False,True,True,True))
# GPy.util.plot.gpplot(Xnew,m+mean, lower+mean, upper+mean)
# pb.plot([X[-1],X[-1]],[-10+mean,15+mean],'k:')
# year_ticks()
# pb.savefig('CO2a-ma.pdf')

# ########################################################################
# ########################################################################
# ## Fig 1b

# ##############################################
# ## create GP model for periodic and aperiodic part

# kp = GPy.kern.periodic_Matern32(D=1,variance=var,lengthscale=lenscl,period=per,n_freq=N,lower=a,upper=b)
# ka1 = GPy.kern.Matern32(D=1,variance=var,lengthscale=lenscl)
# ka2 = GPy.kern.periodic_Matern32(D=1,variance=var,lengthscale=lenscl,period=per,n_freq=N,lower=a,upper=b)
# k = kp + ka1 + ka2 #+ noise

# mm = GPy.models.GP_regression(X,Y,kernel = k)
# mm.tie_params(np.array((3,5)))
# mm.tie_params(np.array((4,6)))

# # mm.constrain_positive(np.array((0,1,3,4,6)))
# # mm.constrain_negative('periodic_Mat32_1_variance')
# # mm.constrain_fixed('_period')
# # mm.constrain_fixed('noise_variance',1e-2)

# #var
# mm.constrain_positive(np.array((0,3)))
# mm.constrain_negative('periodic_Mat32_1_variance')
# #mm.constrain_positive('noise_variance')
# mm.constrain_bounded('noise_variance',0.001,0.2)
# #len
# mm.constrain_bounded(np.array((1,4,6)),0.5,20)
# #per
# mm.constrain_fixed('_period')


# print mm
# mm.checkgrad(verbose=True)
# mm.optimize()
# #mm.optimize_restarts(5)
# print mm

# ##############################################
# # plots

# pb.figure(figsize=(4,4))
# pb.plot(Xnew,Ynew+mean,'r-', label='line 1', linewidth=1.5)
# m, tmp, lower, upper = mm.predict(Xnew)
# GPy.util.plot.gpplot(Xnew,m+mean, lower+mean, upper+mean)
# pb.plot([X[-1],X[-1]],[-10+mean,15+mean],'k:')
# year_ticks()
# pb.savefig('CO2b-t.pdf')

# pb.figure(figsize=(4,4))
# m, tmp, lower, upper = mm.predict(Xnew, slices=(True,False,False,False))
# GPy.util.plot.gpplot(Xnew,m+mean, lower+mean, upper+mean)
# pb.plot([X[-1],X[-1]],[-10+mean,15+mean],'k:')
# year_ticks()
# pb.savefig('CO2b-mp.pdf')

# pb.figure(figsize=(4,4))
# m, tmp, lower, upper = mm.predict(Xnew, slices=(False,True,True,True))
# GPy.util.plot.gpplot(Xnew,m+mean, lower+mean, upper+mean)
# pb.plot([X[-1],X[-1]],[-10+mean,15+mean],'k:')
# year_ticks()
# pb.savefig('CO2b-ma.pdf')


# #

# ## model test
# #mp,vp,inf,sup = mm.predict(Xnew)
# #np.std((mp - Ynew)/np.sqrt(vp))
