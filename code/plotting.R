#################################################
## Libraries and functions used in the code
#################################################
library("R.matlab")
library("latex2exp")
library("foreach")
library("doParallel")
library("parallel")
library("MASS")
numCores <- detectCores()
registerDoParallel(numCores-1)  # use multicore, 8 in total 
library("Rcpp")
library("RcppArmadillo")
library("SQUAREM")
library("beepr")

setwd("C:/Users/au543194/Documents/gitcode/sampleSFS_paper") # working directory
source("code/AllFunctions.R") # Functions to use
###############################################
## Loading data
###############################################
load("data/LACA24.RData")    # data to use
G = ncol(V)
K = nrow(V)
M = V[,order(colSums(V))]
N = 3
################ Compute one solution ##########################
fit = NMFPois(M, N = N, seed = 1:20, tol = 1e-8) 

##########################################################
## Test different signatures and beta simultaneously
##########################################################
maxIter = 10^6
N.sig = c(2:10)           # different number of signatures
each = 50                 # times to run each example
beta = c(0.1,0.5,1)       # different values of beta
beta.vec = rep(beta,each*length(N.sig))         # vector of beta values
sig.vec = rep(N.sig, each = each*length(beta))  # vector of signature number
fit.idx = rep(1:length(N.sig), each = each*length(beta)) # vector of index for fits 
x = c(1:(each*length(N.sig)*length(beta)))               # index of run

# Finds global minimum of PE from five initialization setting a N equal to x
runNMF = function(x){
  print(x)
  return(NMFPois(M, N = x, seed = sample(1:100,5), tol = 1e-5))
}

# Running the function for chosen number of signatures
fits = lapply(N.sig,runNMF)

# Finds the SFS for certain combinations of number of signatures and beta
findsfs = function(x){
  print(x)
  start_time = Sys.time()
  output = sampleSFS(P = fits[[fit.idx[x]]]$P, E = fits[[fit.idx[x]]]$E, maxIter = maxIter, beta = beta.vec[x], check = 1000, eps = 1e-8)
  end_time = Sys.time()
  output$time.elap = end_time - start_time
  return(output)
}
## Running the function 'each' number of times for each combination
system.time({
resultBRCA = lapply(x,findsfs)
})

## creates a data.frame of the results with summary of the observations for all combinations
varSummary = function(variable, list, variable.descrip, data.descrip){
  size = sapply(x, function(x) list[[x]][[variable]]) 
  mean = tapply(size, paste(sig.vec,beta.vec), mean)           # mean of size for each combination
  min = tapply(size, paste(sig.vec,beta.vec), min)             # min of size for each combination
  max = tapply(size, paste(sig.vec,beta.vec), max)             # max of size for each combination
  sd = tapply(size, paste(sig.vec,beta.vec), sd)               # standard deviation of size for each combination
  q25 = tapply(size, paste(sig.vec,beta.vec), function(x) quantile(x, probs = 0.25))
  q75 = tapply(size, paste(sig.vec,beta.vec), function(x) quantile(x, probs = 0.75))
  
  ## data frame of average change of the entries in P
  data = data.frame(beta.val = rep(paste("beta = ",beta),length(N.sig)), sig.no = rep(c(10,N.sig[-9]), each = length(beta)), 
                        mean = mean, min = min, max = max, sd = sd, id = variable.descrip, Cancer = data.descrip,
                    q25 = q25, q75 = q75)
  data$beta.val = factor(data$beta.val)
  return(data)
}
load("results/resultBRCA.RData")
load("results/resultLACA.RData")
# get data about size of SFS
dataBRCAchange = varSummary(variable = "avgChangeFinal", 
                            list = resultBRCA, variable.descrip = "Size of SFS", 
                            data.descrip = "Breast Cancer")
dataLACAchange = varSummary(variable = "avgChangeFinal", 
                            list = resultLACA, variable.descrip = "Size of SFS", 
                            data.descrip = "Lung A. Cancer")

data3 = rbind(dataBRCAchange,dataLACAchange)
data3$id = factor(data3$id)
data3$Cancer = factor(data3$Cancer)

## Plot of the results 
ggplot(data3, aes(x = sig.no + rep(c(-0.06,0,0.06),18), y = mean, group = beta.val, color = beta.val))+
  geom_point()+geom_line()+
  facet_wrap(~ Cancer, scales = "free")+
  geom_errorbar(aes(ymin = q25, ymax = q75), size = 0.6, width = 0.2)+
  ylab(expression(italic(avg*symbol("\341")*P^S*symbol("\361"))))+
  xlab("N (Number of mutational processes)")+theme_bw()+
  labs(color = " ", group = " ", fill = " ")+
  theme(legend.position = "bottom", text = element_text(size=15))+
  scale_x_continuous(breaks = c(2,3,4,5,6,7,8,9,10))

# number of iterations and time
dataLACAiter = varSummary(variable = "totalIter", list = resultLACA[beta.vec == 0.5], 
                          variable.descrip = "Iterations before stopping", 
                          data.descrip = "Lung A. Cancer")
dataLACAtime = varSummary(variable = "time.elap", list = resultLACA[beta.vec == 0.5], 
                          variable.descrip = "Time before stopping", 
                          data.descrip = "Lung A. Cancer")
dataBRCAiter = varSummary(variable = "totalIter", list = resultBRCA[beta.vec == 0.5], 
                          variable.descrip = "Iterations before stopping", 
                          data.descrip = "Breast Cancer")
dataBRCAtime = varSummary(variable = "time.elap", list = resultBRCA[beta.vec == 0.5], 
                          variable.descrip = "Time before stopping", 
                          data.descrip = "Breast Cancer")

data3 = rbind(dataBRCAiter,dataLACAiter)
data3$id = factor(data3$id)
data3$Cancer = factor(data3$Cancer)

t1 = ggplot(data3, aes(x = sig.no, y = mean, col = Cancer))+
  geom_point()+geom_line()+
  #facet_grid(cols = vars(Cancer), scales = "free_y" ,shrink = TRUE)+
  geom_errorbar(aes(ymin = q25, ymax = q75), size = 0.6, width = 0.2)+
  ylab("Iterations before stopping")+xlab("Number of mutational processes")+theme_bw()+
  labs(color = " ", group = " ", fill = " ")+
  theme(legend.position = "bottom", text = element_text(size=10))+
  scale_x_continuous(breaks = c(2,3,4,5,6,7,8,9,10))
t1 

data3 = rbind(dataBRCAtime,dataLACAtime)
data3$id = factor(data3$id)
data3$Cancer = factor(data3$Cancer)

t2 = ggplot(data3, aes(x = sig.no, y = mean, col = Cancer))+
  geom_point()+geom_line()+
  #facet_grid(cols = vars(Cancer), scales = "free_y" ,shrink = TRUE)+
  geom_errorbar(aes(ymin = q25, ymax = q75), size = 0.6, width = 0.2)+
  ylab("Time before stopping(sec)")+xlab("Number of mutational processes")+theme_bw()+
  labs(color = " ", group = " ", fill = " ")+
  theme(legend.position = "bottom", text = element_text(size=10))+
  scale_x_continuous(breaks = c(2,3,4,5,6,7,8,9,10))
t2

ggarrange(t1, t2, common.legend = T, legend = "bottom")
plot(c(10,2:9),dataBRCAtime$mean)
points(c(10,2:9),dataBRCAiter$mean*c(10,2:9)/100000)
#################################################################
## plot to compare with SVD area from FACPACK and change of beta
#################################################################
## result from SFS
D = fit$P%*%fit$E
# writeMat(con = "BRCAN3.mat",D = D) # export to use on FAC-PACK

## results from sampling algorithm
sample = 15000
sfsres = sampleSFS(P = fit$P, E = fit$E, maxIter = sample, check = sample, beta = 0.5, eps = 1e-8)
res = samplesToSVD(Presults = sfsres$P_lastCheckResults, Eresults = sfsres$E_lastCheckResults, N = 3)

#### Results from polygon inflation algorithm 
facpack.res = readMat(con = "data/results_facpack_all_LACAN3_transpose.mat") # Choose the right results from FAC-PACK
AFS = rbind(facpack.res$AFS[[1]][[1]],facpack.res$AFS[[3]][[1]],facpack.res$AFS[[2]][[1]])
type = c(rep(1,nrow(facpack.res$AFS[[1]][[1]])),rep(2,nrow(facpack.res$AFS[[3]][[1]])),rep(3,nrow(facpack.res$AFS[[2]][[1]])))

dat.T = data.frame(x = res$E.points[,1], y = -res$E.points[,2])  # Choose PE.points or ET.points
dat.Pol = data.frame(x = AFS[,1], y = AFS[,2], type = factor(type))

#### PLOT 
ggplot(dat.T, aes(x = x, y = y))+
  geom_point(size = 0.2, col = c(rep("darkgrey", nrow(dat.T)-500), rep("black",500)), shape = 3)+
  geom_polygon(data = dat.Pol, aes(x = x, y = y, group = type, colour = type), 
                
               fill=NA, size = 1.2, linetype = "twodash")+
  labs(x = expression(alpha[1]), y = expression(alpha[2]))+
  theme_bw()+
  theme(text = element_text(size = 15), legend.position = "none")+
  geom_segment(aes(x = 0.8, y = 0.3, xend = 0.05, yend = 0.02),
               arrow = arrow(length = unit(0.3, "cm")), colour = "#E69F00", size = 1.1)+
  scale_color_manual(values = c("#E69F00","#56B4E9","#009E73"), name = element_blank(), labels = c("Signature 1", "Signature 2", "Signature 3") )
## only plot 1000 points of signature 2 for different beta values
plots = list()
sample = 1000
beta = c(0.1,0.5,1)
for(b in 1:length(beta)){
  sfsres = sampleSFS(P = fit$P, E = fit$E, maxIter = sample, check = sample, beta = beta[b], eps = 1e-8)
  resT = samplesToSVD(Presults = sfsres$P_lastCheckResults, Eresults = sfsres$E_lastCheckResults, N = 3)
  dat.T = data.frame(x = resT$P.points[,1], y = -resT$P.points[,2])
  alpha.val = seq(0,1, length.out = sample)
  
  plots[[b]] = ggplot(dat.T, aes(x = x, y = y))+
    geom_point(size = 1, alpha = 0.6)+
    labs(x = expression(alpha[1]), y = expression(alpha[2]))+
    theme_bw()+
    lims(x = c(0.3,1), y = c(0.8,1.6))+
  theme(legend.position = "none", text = element_text(size = 15))
}

plots[[3]]
#########################################################################
#### Signature plot of SFS area from sampling algorithm
##########################################################################
dev1 = sampleSFS(P=fit$P, E = fit$E, maxIter = 10^5, beta = 0.5, check = 1000)
prob.min = dev1$Pminimum
prob.max = dev1$Pmaximum

# Converting into right data frame
sig = rep(c(1:N), dev1$totalIter)
mut <- c("C > A","C > G","C > T","T > A","T > C","T > G")
sub = rep(mut, each = 16)
dat1 = data.frame(m = factor(rownames(V), levels=unique(rownames(V))), sub = sub, S = t(prob.min))
datmin = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long", v.names = "min")
dat1 = data.frame(m = factor(rownames(V), levels=unique(rownames(V))),S = t(prob.max))
datmax = reshape(dat1, varying = colnames(dat1)[-1], direction = "long", v.names = "max")
data2 = merge(datmin,datmax, by = c("m","time","id"))
data2$time = factor(data2$time)

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}
## Plotting the signatures
#col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")
g1 = ggplot(data2, aes(x = m, y = min))+
  #geom_errorbar(aes(ymin = min, ymax = min, col = time), width = 1.05, size = 1)+
  #geom_errorbar(aes(ymin = max, ymax = max, col = time), width = 1.05, size = 1)+
  geom_bar(aes(x = m, y = max), stat = "identity", width = 0.78, fill = "tomato2")+
  #geom_errorbar(aes(ymin = min, ymax = max, col = time), lwd = 1, size = 0.5, linetype = "11")+
  geom_bar(stat = "identity", width = 0.8, fill = "grey")+
  facet_grid(rows = vars(time), cols = vars(sub), scales = "free", switch = "x")+theme_bw()+
  theme(text = element_text(size=12, face = "bold"), axis.text.x=element_blank(),axis.text.y = element_text(size = 8),axis.ticks = element_blank(), 
        legend.position = "none",strip.text.y.right = element_text(angle = 0, size = 15), panel.spacing.x = unit(0,"line"),
        strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
        strip.text.x = element_text(size = 9))+ 
  ylab("Probability")+xlab("Mutation types")+ggtitle("Signatures")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=equal_breaks(n=3, s=0.2), 
                     expand = c(0.05, 0))
  #scale_colour_grey()+
  #scale_fill_grey()
g1

# include colors on strips for N = 3
g11 <- ggplot_gtable(ggplot_build(g1))
stripr <- which(grepl('strip-r', g11$layout$name))
fills <- c("#E69F00", "#56B4E9", "#009E73")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g11$grobs[[i]]$grobs[[1]]$childrenOrder))
  g11$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

####################################################################
## exposure plot
##################################################################

## equivalent to signature plot
prob.min = dev1$Eminimum%*%diag(1/colSums(dev1$E_lastCheckResults[1:N,]))
prob.max = dev1$Emaximum%*%diag(1/colSums(dev1$E_lastCheckResults[1:N,]))
dat1 = data.frame(m = factor(colSums(M), levels=unique(colSums(M))), S = t(prob.min))
datmin = reshape(dat1, varying = colnames(dat1)[-1], direction = "long", v.names = "min")
dat1 = data.frame(m = factor(colSums(M), levels=unique(colSums(M))),S = t(prob.max))
datmax = reshape(dat1, varying = colnames(dat1)[-1], direction = "long", v.names = "max")

dat2 = merge(datmin,datmax, by = c("m","time","id"))
dat2$time = factor(dat2$time)

## different plots 
g2 = ggplot(dat2, aes(x = m, y = min))+
  #geom_errorbar(aes(ymin = min, ymax = min, col = time), width = 1, size = 1)+
  #geom_errorbar(aes(ymin = max, ymax = max, col = time), width = 1, size = 1)+
  geom_bar(aes(x = m, y = max), stat = "identity", width = 0.78, fill = "tomato2")+
  #geom_errorbar(aes(ymin = min, ymax = max, col = time), lwd = 1, size = 0.5, linetype = "11")+
  geom_bar(stat = "identity", width = 0.8, fill = "grey")+
  facet_grid(cols = vars(time))+
  theme_bw()+
  theme(text = element_text(size=12, face = "bold"), axis.text.x=element_text(angle = 315, size = 5, hjust = 0.7, vjust = -0.6), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "none",strip.text.y.right = element_text(angle = 0))+ 
  ylab("Probability")+xlab("Patients")+ggtitle("   Normalized Exposures")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0,0.5,1))+
  coord_flip()
  #scale_colour_grey()+
  #scale_fill_grey()

g2

# include colors on strips for N = 3
g22 <- ggplot_gtable(ggplot_build(g2))
stripr <- which(grepl('strip-t', g22$layout$name))
fills <- c("#E69F00", "#56B4E9", "#009E73")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g22$grobs[[i]]$grobs[[1]]$childrenOrder))
  g22$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# final plot
ggarrange(g11,g22, labels = c("(A)","(B)"), widths = c(2, 1))

########################################################################
## Initialization 
########################################################################
each = 100
noInit = 10
# running the algorithm 1000 times with different initializations
NMFres = foreach(i = 1:(each*noInit)) %dopar% {
  library(SQUAREM)
  NMFPois(M,N = 3, seed = i+2000, tol = 1e-6)
}

## plot of increase in initialization
GKLmat = matrix(sapply(1:(noInit*each), function(x) NMFres[[x]]$gkl), nrow = noInit, ncol = each)
GKLmatmin = apply(GKLmat,2,cummin)
GKLavg = apply(GKLmatmin,1,mean)
GKLsd = apply(GKLmatmin,1,sd)
GKL0.05 = apply(GKLmatmin,1, function(x) quantile(x, probs = 0.05))
GKL0.95 = apply(GKLmatmin,1, function(x) quantile(x, probs = 0.95))


gkldat = data.frame(avg = GKLavg, sd = GKLsd, q05 = GKL0.05, 
                  q95 = GKL0.95, 
                  Init = c(1:noInit), Cancer = "Breast Cancer") # Cancer can also be set to "Lung A. Cancer"

ggplot(gkldat, aes(x = Init, y = avg, col = Cancer))+
  geom_point()+geom_line()+
  geom_errorbar(aes(ymin = q05, ymax = q95), size = 0.6, width = 0.2)+
  ylab("minimum D(M|PE)")+xlab("Number of initializations")+theme_bw()+
  labs(color = " ", group = " ", fill = " ")+
  theme(legend.position = "none", text = element_text(size = 12))+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))


## plot of the change in the SVD area with three initializations
minrun = apply(GKLmat, 2, which.min) # finding which of the ten initializations gave the minimum 
minimum = apply(GKLmat,2, min)
# extracting best 100 fits
Pmat = sapply(1:each, function(x) NMFres[[noInit*(x-1) + minrun[x]]]$P)
Pmat = matrix(Pmat, nrow = 96)
Emat = sapply(1:each, function(x) t(NMFres[[noInit*(x-1) + minrun[x]]]$E))
Emat = matrix(Emat, nrow = 24)

res = samplesToSVD(Presults = t(Pmat), Eresults = t(Emat), N = 3)

#### Results from polygon inflation algorithm 
facpack.res = readMat(con = "data/results_facpack_all_BRCAN3.mat") # Choose the right results from FAC-PACK
AFS = rbind(facpack.res$AFS[[1]][[1]],facpack.res$AFS[[3]][[1]],facpack.res$AFS[[2]][[1]])
type = c(rep(1,nrow(facpack.res$AFS[[1]][[1]])),rep(2,nrow(facpack.res$AFS[[3]][[1]])),rep(3,nrow(facpack.res$AFS[[2]][[1]])))

dat.T = data.frame(x = res$P.points[,1], y = -res$P.points[,2])  # Choose PE.points or ET.points
dat.Pol = data.frame(x = AFS[,1], y = AFS[,2], type = factor(type))

#### PLOT 
ggplot(dat.T, aes(x = x, y = y))+
  geom_polygon(data = dat.Pol, aes(x = x, y = y, group = type), fill="grey", size = 1.2)+
  geom_point(size = 1, shape = 3)+
  labs(x = expression(alpha[1]), y = expression(alpha[2]))+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 15))

## check for same optimal factorization
estimates = lapply(idx.min, function(x) NMFres[[noInit*(x-1) + minrun[x]]]$P%*%NMFres[[noInit*(x-1) + minrun[x]]]$E)
dist = numeric(length(estimates))
for(i in 1:length(estimates)){
  dist[i] = sum((D - estimates[[i]])^2)
}


############################################################
### Parametric bootstrap results
###########################################################
iter = 100
res1 = boot_pois(fit$P,fit$E, iter, same.init = FALSE)
res2 = boot_pois(fit$P,fit$E, iter, same.init = TRUE)

## Converted to SVD area 
resT1 = samplesToSVD(Presults = res1$Presults, Eresults = res1$Eresults, N = 3, Mfit = fit$P%*%fit$E)
resT2 = samplesToSVD(Presults = res2$Presults, Eresults = res2$Eresults, N = 3, Mfit = fit$P%*%fit$E)

dat.T = data.frame(x = resT2$P.points[,1], y = resT2$P.points[,2]) # change to either resT1 or resT2

ggplot(dat.T, aes(x = x, y = y))+
  geom_polygon(data = dat.Pol, aes(x = x, y = y, group = type), fill="grey", size = 1.2)+
  geom_point(size = 1, shape = 3)+
  labs(x = expression(alpha[1]), y = expression(alpha[2]))+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 15))
