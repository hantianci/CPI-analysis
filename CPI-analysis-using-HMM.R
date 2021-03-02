library(astsa)
library(tseries)
library(depmixS4)
library(TSA)

# Load Data
cpi = structure(c(216.687, 216.741, 217.631, 218.009, 218.178, 217.965,
                  218.011, 218.312, 218.439, 218.711, 218.803, 219.179,
                  220.223, 221.309, 223.467, 224.906, 225.964, 225.722,
                  225.922, 226.545, 226.889, 226.421, 226.230, 225.672,
                  226.665, 227.663, 229.392, 230.085, 229.815, 229.478,
                  229.104, 230.379, 231.407, 231.317, 230.221, 229.601,
                  230.280, 232.166, 232.773, 232.531, 232.945, 233.504,
                  233.596, 233.877, 234.149, 233.546, 233.069, 233.049,
                  233.916, 234.781, 236.293, 237.072, 237.900, 238.343,
                  238.250, 237.852, 238.031, 237.433, 236.151, 234.812,
                  233.707, 234.722, 236.119, 236.599, 237.805, 238.638,
                  238.654, 238.316, 237.945, 237.838, 237.336, 236.525,
                  236.916, 237.111, 238.132, 239.261, 240.229, 241.018,
                  240.628, 240.849, 241.428, 241.729, 241.353, 241.432,
                  242.839, 243.603, 243.801, 244.524, 244.733, 244.955,
                  244.786, 245.519, 246.819, 246.663, 246.669, 246.524,
                  247.867, 248.991, 249.554, 250.546, 251.588, 251.989,
                  252.006, 252.146, 252.439, 252.885, 252.038, 251.233,
                  251.712, 252.776, 254.202, 255.548, 256.092, 256.143,
                  256.571, 256.558, 256.759, 257.346, 257.208, 256.974,
                  257.971, 258.678, 258.115, 256.389, 256.394, 257.797,
                  259.101, 259.918, 260.280),
                .Tsp = c(2010, 2020+8/12, 12), class = "ts")
dput(cpi)
par(mfrow=c(1,1))
plot(cpi)
# Box-Cox
BoxCox.ar(cpi, method='burg')
# Check Stationarity
acf(cpi)
pacf(cpi)
adf.test(cpi)
pp.test(cpi) 
acf(diff(cpi))
pacf(diff(cpi))
adf.test(diff(cpi))
pp.test(diff(cpi)) 
plot(diff(cpi))
# Hidden Markov Model (HMM)
d = diff(cpi)
y = ts(d, start=2010+1/12, freq=12) # make data depmix friendly
mod3 <- depmix(y~1, nstates=3, data=data.frame(y))
set.seed(100)
summary(fm3 <- fit(mod3))
# Graphics of CPI
plot(cpi, main="", ylab='CPI', col=gray(.7),
     ylim=c(210,265))
culer = 4-posterior(fm3)[,1] # switch labels 1 and 3
text(cpi, col=culer, labels=4-posterior(fm3)[,1])
# Graphics of diff CPI
layout(matrix(c(1,2, 1,3), 2), heights=c(1,.75))
par(mar=c(2.5,2.5,.5,.5), mgp=c(1.6,.6,0))
plot(y, main="", ylab='First Order Difference CPI', col=gray(.7),
     ylim=c(-2,2.3))
culer = 4-posterior(fm3)[,1] # switch labels 1 and 3
text(y, col=culer, labels=4-posterior(fm3)[,1])
# MLEs
para.mle = as.vector(getpars(fm3)[-(1:3)])
permu = matrix(c(0,0,1,0,1,0,1,0,0), 3,3) # for the label switch
(mtrans.mle = permu%*%round(t(matrix(para.mle[1:9],3,3)),3)%*%permu)
(norms.mle = round(matrix(para.mle[10:15],2,3),3)%*%permu)
acf(y^2, panel.first=grid(lty=2) )
hist(y, 25, prob=TRUE, main='')
culer=c(1,2,3); pi.hat = colSums(posterior(fm3)[-1,2:4])/length(y)
for (i in 1:3) { mu=norms.mle[1,i]; sig = norms.mle[2,i]
x = seq(-1.8,2.2, by=.001)
lines(x, pi.hat[4-i]*dnorm(x, mean=mu, sd=sig), col=culer[i]) }
# Bootstrap
set.seed(666); n.obs = length(y); n.boot = 20
para.star = matrix(NA, nrow=n.boot, ncol = 15)
respst <- para.mle[10:15]; trst <- para.mle[1:9]
for ( nb in 1:n.boot ){
  mod <- simulate(mod3)
  y.star = as.vector(mod@response[[1]][[1]]@y)
  dfy = data.frame(y.star)
  mod.star <- depmix(y.star~1, data=dfy, respst=respst, trst=trst, nst=3)
  fm.star = fit(mod.star, emcontrol=em.control(tol = 1e-5), verbose=FALSE)
  para.star[nb,] = as.vector(getpars(fm.star)[-(1:3)]) }
# bootstrap stnd errors
SE = sqrt(apply(para.star,2,var) + (apply(para.star,2,mean)-para.mle)^2)
(SE.mtrans.mle = permu%*%round(t(matrix(SE[1:9],3,3)),3)%*%permu)
(SE.norms.mle = round(matrix(SE[10:15], 2,3),3)%*%permu)

