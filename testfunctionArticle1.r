library(spatstat)
library(lattice)
library(sp)
library(maptools)
library(raster)
library(geostatsp)
library(rgdal)
library(lattice)
library(caret)
library(rgeos)
library(scales)


#------------------------------------------------------------------------------
# 								Simulation data
#------------------------------------------------------------------------------

# Set up some data
# 1 #  Set up data.ppp, cov.list, ppmform and quads
# Generate XY grid
set.seed(10013)
XY = expand.grid(seq(0, 100, 1), seq(0, 100, 1))
X = XY[,1]
Y = XY[,2]

# Generate 2 covariates for PPM

v1 = (X - 30)^2 + (Y - 70)^2 - 0.5*X*Y

v2 = (X - 70)^2 + (Y - 60)^2 + 0.9*X*Y 


#levelplot(v1 ~ X + Y)
#levelplot(v2 ~ X + Y)


v1 = -1*scale(v1)
v2 = -1*scale(v2)

# Matrix of covariates
vmat = as.matrix(data.frame(1, v1, v1^2, v2, v2^2))

# Generate true PPM coefficients based on linear and quadratic terms for 2 covariates and including bias
sp1_coef = c(-6.5, 4, -1, 2, -0.6)
sp1_int = exp(vmat %*% sp1_coef)

sp2_coef = c(-4.4, 1.8, -1, 1.5, -0.9)
sp2_int = exp(vmat %*% sp2_coef)

sp3_coef = c(-3.5, -0.5, -0.8, 1, -0.8)
sp3_int = exp(vmat %*% sp3_coef)


sp_int.list = list(sp1_int, sp2_int, sp3_int)
  
# Plot the intensities created

levelplot(sp1_int ~ X + Y)
levelplot(sp2_int ~ X + Y)
levelplot(sp3_int ~ X + Y)

# Create pixel images of intensity surfaces for spatstat

sp1_int_im = as.im(data.frame(x = X, y = Y, z = sp1_int))
sp2_int_im = as.im(data.frame(x = X, y = Y, z = sp2_int))
sp3_int_im = as.im(data.frame(x = X, y = Y, z = sp3_int))


# Simulate species patterns

sp1_sim = rpoispp(sp1_int_im)
sp2_sim = rpoispp(sp2_int_im)
sp3_sim = rpoispp(sp3_int_im)

sp1_sim
sp2_sim
sp3_sim

plot(sp1_sim, cex = 0.6)
plot(sp2_sim, add = TRUE, col = "red", cex = 0.6)
plot(sp3_sim, add = TRUE, col = "blue", cex = 0.6)

sp_sim.list = list(sp1_sim, sp2_sim, sp3_sim)

# Look at the correlation between intensity surfaces
#all
cor1_2 = cor(as.vector(sp1_int), as.vector(sp2_int), use = "complete.obs")
cor1_3 = cor(as.vector(sp1_int), as.vector(sp3_int), use = "complete.obs")
cor2_3 = cor(as.vector(sp2_int), as.vector(sp3_int), use = "complete.obs")


## Create list of coavriates

cov.list = list()
for (v in 1:4)
{
v.v = as.im(data.frame(x = X, y = Y, z = vmat[,(v + 1)]))
cov.list[[v]] = v.v
}
names(cov.list) = c("v1", "v1.2", "v2", "v2.2")

# set up model formula
cov.mat = vmat[,2:5]
ppmform = as.formula(paste("~", paste(colnames(cov.mat), collapse = "+")))


##################################################
# 				Test new simulation 
#-------------------------------------------------

# Call the different functions needed for the test
source("functionTestsimarticle1.R")

#_____________________________________________________________________________

# Create a confusion matrix from the given outcomes, whose rows correspond
# to the actual and the columns to the predicated classes.
createConfusionMatrix <- function(act, pred) {
pred <- pred[order(act)]
act  <- act[order(act)]
sapply(split(pred, act), tabulate, nbins=3)
}


# Function to combine the different methods to compare

Testsims = function(hidepct, n.sims, sp_sim.list, sp_int_im, n.sp=n.sp, k = k, ks=ks, nstart=nstart, cov.list, cov.bias=NULL, kVal=NULL, kAreaInt=NULL, delta_max=delta_max, delta_min=delta_min, delta_step =delta_step, num.add = num.add)
{
  quads = ppp(X, Y, window = win)
  
  RSSknn = meanRSSknn = IMSEknn = RSSkmeans = meanRSSkmeans = IMSEkmeans =  RSSkps = meanRSSkps = IMSEkps = 
    RSSrand = meanRSSrand = IMSErand = RSSequal = meanRSSequal = IMSEequal =  RSSCF = meanRSSCF = IMSECF = 
    RSSindiv = meanRSSindiv = IMSEindiv =RSSLoopT = meanRSSLoopT = IMSELoopT =
    RSSLoopE = meanRSSLoopE = IMSELoopE = RSSLoopA = meanRSSLoopA = IMSELoopA=
    sumcorknn1 = sumcorkmeans1 = sumcorrand1 = sumcorequal1 = sumcorindiv1 = sumcorLoopE1 =
    sumcorLoopA1 = sumcorLoopT1 = sumcorkps1 = sumcorCF1 = 
    sumcorknn2 = sumcorkmeans2 = sumcorrand2 = sumcorequal2 = sumcorindiv2 = sumcorLoopE2 =
    sumcorLoopA2 = sumcorLoopT2 = sumcorknn3 = sumcorkmeans3 = sumcorrand3 = sumcorequal3 = sumcorindiv3 = sumcorLoopE3 =
    sumcorLoopA3 = sumcorLoopT3 = sumcorkps2 = sumcorCF2 = sumcorkps3 = sumcorCF3 = matrix(NA, n.sims, length(hidepct))
  accmatknn = accmatkmeans = accmatrand = accmatequal = accmatindiv = accmatLoop = 
    accmatLoopA = accmatLoopT = accmatLoopE = accmatkps = accmatCF = matrix(NA, n.sims, length(hidepct))
  
  
  knnpred = kmeanspred = randpred = equalpred = LoopApred = LoopTpred = LoopEpred =
    indivpred = kpspred = CFpred = array(NA, c(quads$n, 3, n.sims, length(hidepct)))
  
  coef.knn.mat = coef.kmeans.mat = coef.rand.mat = coef.eq.mat = coef.kps.mat = coef.CF.mat = 
    array(NA, c(15, 1, n.sims, length(hidepct)))
  
  coef.LA.mat = coef.LT.mat = coef.LE.mat =  coef.ind.mat = 
    array(NA, c(15, 1, n.sims, length(hidepct)))
  
  W.knn = W.kmeans = W.rand = W.kps = W.CF = W.ind = W.LA = W.LT = W.LE = replicate(3, list())
  sp.true = replicate(3, list())
  
  if(length(hidepct)==1){
    se.knn1 = se.kmeans1 = se.rand1 = se.CF1 = se.kps1 = se.LA1 = se.LT1 = se.LE1 = se.indiv1 = list()
  }else{
    se.knn1 = se.kmeans1 = se.rand1 = se.CF1 = se.kps1 = se.LA1 = se.LT1 = se.LE1 = se.indiv1 = rep(list(list()), length(hidepct))
  }
  
  if(length(hidepct)==1){
    sew.knn1 = sew.kmeans1 = sew.rand1 = sew.CF1 = sew.kps1 = sew.LA1 = sew.LT1 = sew.LE1 = sew.indiv1 = list()
  }else{
    sew.knn1 = sew.kmeans1 = sew.rand1 = sew.CF1 = sew.kps1 = sew.LA1 = sew.LT1 = sew.LE1 = sew.indiv1 = rep(list(list()), length(hidepct))
  }
  
  if(length(hidepct)==1){
    se.knn2 = se.kmeans2 = se.rand2 = se.CF2 = se.kps2 = se.LA2 = se.LT2 = se.LE2 = se.indiv2 = list()
  }else{
    se.knn2 = se.kmeans2 = se.rand2 = se.CF2 = se.kps2 = se.LA2 = se.LT2 = se.LE2 = se.indiv2 = rep(list(list()), length(hidepct))
  }
  
  if(length(hidepct)==1){
    sew.knn2 = sew.kmeans2 = sew.rand2 = sew.CF2 = sew.kps2 = sew.LA2 = sew.LT2 = sew.LE2 = sew.indiv2 = list()
  }else{
    sew.knn2 = sew.kmeans2 = sew.rand2 = sew.CF2 = sew.kps2 = sew.LA2 = sew.LT2 = sew.LE2 = sew.indiv2 = rep(list(list()), length(hidepct))
  }
  
  if(length(hidepct)==1){
    se.knn3 = se.kmeans3 = se.rand3 = se.CF3 = se.kps3 = se.LA3 = se.LT3 = se.LE3 = se.indiv3 = list()
  }else{
    se.knn3 = se.kmeans3 = se.rand3 = se.CF3 = se.kps3 = se.LA3 = se.LT3 = se.LE3 = se.indiv3 = rep(list(list()), length(hidepct))
  }
  
  if(length(hidepct)==1){
    sew.knn3 = sew.kmeans3 = sew.rand3 = sew.CF3 = sew.kps3 = sew.LA3 = sew.LT3 = sew.LE3 = sew.indiv3 = list()
  }else{
    sew.knn3 = sew.kmeans3 = sew.rand3 = sew.CF3 = sew.kps3 = sew.LA3 = sew.LT3 = sew.LE3 = sew.indiv3 = rep(list(list()), length(hidepct))
  }
  
  for (i in 1:length(hidepct))
  {
    pct_hidden = hidepct[i]
    
    for (j in 1:n.sims)
    {
      # hide some observations
      sp_hide.list = sp_sub.list = train.list = sp_test.list = list()
      coordtestx.list = coordtesty.list = markshide.list = markstest.list = list()
      coordsubx.list = coordsuby.list = marksub.list = list()
      
      for (l in 1:n.sp) {
        sp_hide.list[[l]] = sample(1:sp_sim.list[[l]]$n, floor(pct_hidden*sp_sim.list[[l]]$n))
        sp_sub.list[[l]]  = sp_sim.list[[l]][-sp_hide.list[[l]]]
        train.list[[l]]   = ppp(x = sp_sub.list[[l]]$x, y = sp_sub.list[[l]]$y, window = win)
        sp_test.list[[l]] = sp_sim.list[[l]][sp_hide.list[[l]]]
        
        coordtestx.list[[l]] = sp_test.list[[l]]$x
        coordtesty.list[[l]] = sp_test.list[[l]]$y
        markshide.list[[l]] = rep(paste("Hidden", l, sep = ""), sp_test.list[[l]]$n)
        markstest.list[[l]] = rep(paste("Sp", l, sep = ""), sp_test.list[[l]]$n)
        
        coordsubx.list[[l]] = sp_sub.list[[l]]$x
        coordsuby.list[[l]] = sp_sub.list[[l]]$y
        marksub.list[[l]] = rep(paste("Sp", l, sep = ""), sp_sub.list[[l]]$n)
        
        l=l+1
      }
      sp.true[[i]][[j]] = sp_hide.list
      
      all_test = ppp(x = c(unlist(coordtestx.list)), 
                     y = c(unlist(coordtesty.list)), window = win,
                     marks = c(unlist(markshide.list)))
      
      all_test2 = ppp(x = c(unlist(coordtestx.list)), 
                      y = c(unlist(coordtesty.list)), window = win,
                      marks = c(rep("Unknown", all_test$n)))
      
      test_labels = as.vector(unlist(markstest.list))
      
      all_true = ppp(x = c(unlist(coordsubx.list)), 
                     y = c(unlist(coordsuby.list)), window = win,
                     marks = c(unlist(marksub.list)))
      
      
      datappp = superimpose.ppp(all_true, all_test2)
      
      if(is.null(cov.bias)){
        cov.list = cov.list
      }else{#--- Set observer bias variables to kVal 
        pred.list = cov.list
        set.Val = cov.bias #Variables to set to a certain value
        for (v in set.Val){
          pred.list[[v]]$v = kVal*pred.list[[v]]$v
        }
      }
      
      ###
      #  Mixture model
      ###---
      simknn = ppmMixEngine(datappp = datappp, quads = quads.win, all_true=all_true,
                            all_test=all_test, initweights = "knn", n.sp=n.sp, sp_int_im,
                            k=k, ks=ks, nstart=nstart, ppmform = ppmform, cov.list. = cov.list,
                            cov.bias = cov.bias, kVal = kVal, kAreaInt = kAreaInt,
                            verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simkmeans = ppmMixEngine(datappp = datappp, quads = quads.win, all_true=all_true,
                               all_test=all_test, initweights = "kmeans", n.sp=n.sp,sp_int_im,
                               k=k, ks=ks, nstart=nstart, ppmform = ppmform, cov.list. = cov.list,
                               cov.bias = cov.bias, kVal = kVal, kAreaInt = kAreaInt,
                               verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simrandom = ppmMixEngine(datappp = datappp, quads = quads.win, all_true=all_true,
                               all_test=all_test, initweights = "random", n.sp=n.sp,sp_int_im,
                               k=k, ks=ks, nstart=nstart, ppmform = ppmform, cov.list. = cov.list,
                               cov.bias = cov.bias, kVal = kVal, kAreaInt = kAreaInt,
                               verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simCF = ppmMixEngine(datappp = datappp, quads = quads.win, all_true=all_true, 
                           all_test=all_test, initweights = "CoinF", n.sp=n.sp,sp_int_im,
                           k=NULL, ks=NULL, nstart=NULL, ppmform = ppmform, cov.list. = cov.list,
                           cov.bias = cov.bias, kVal = kVal, kAreaInt = kAreaInt,
                           verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simkps = ppmMixEngine(datappp = datappp, quads = quads.win, all_true=all_true,
                            all_test=all_test, initweights = "kps", n.sp=n.sp,sp_int_im,
                            k=NULL, ks=ks, nstart=nstart, ppmform = ppmform, cov.list. = cov.list,
                            cov.bias = cov.bias, kVal = kVal, kAreaInt = kAreaInt,
                            verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      # for performance measures
      knn_weights = simknn$New_weights[((all_true$n)+1):nrow(simknn$New_weights),]
      pred.knn    = knn_weights
      W.knn[[i]][[j]] = knn_weights
      
      kmeans_weights = simkmeans$New_weights[((all_true$n)+1):nrow(simkmeans$New_weights),]
      pred.kmeans    = kmeans_weights
      W.kmeans[[i]][[j]] = kmeans_weights
      
      random_weights = simrandom$New_weights[((all_true$n)+1):nrow(simrandom$New_weights),]
      pred.random    = random_weights
      W.rand[[i]][[j]] = random_weights
      
      kps_weights = simkps$New_weights[((all_true$n)+1):nrow(simkps$New_weights),]
      pred.kps    = kps_weights
      W.kps[[i]][[j]] = kps_weights
      
      #equal_weights = simequal$New_weights[((all_true$n)+1):nrow(simequal$New_weights),]
      #pred.equal    = equal_weights
      
      CF_weights = simCF$New_weights[((all_true$n)+1):nrow(simCF$New_weights),]
      pred.CF    = CF_weights
      W.CF[[i]][[j]] = CF_weights

      
      accmatknn[j, i] = Accuracy(all_true, knn_weights, test_labels, n.sp)
      RSSknn[j, i] = RSS(pred.knn, test_labels)
      meanRSSknn[j, i] = RSS(pred.knn, test_labels)/length(test_labels)
      
      accmatkmeans[j, i] = Accuracy(all_true, kmeans_weights, test_labels, n.sp)
      RSSkmeans[j, i] = RSS(pred.kmeans, test_labels)
      meanRSSkmeans[j, i] = RSS(pred.kmeans, test_labels)/length(test_labels)
      
      accmatrand[j, i] = Accuracy(all_true, random_weights, test_labels, n.sp)
      RSSrand[j, i] = RSS(pred.random, test_labels)
      meanRSSrand[j, i] = RSS(pred.random, test_labels)/length(test_labels)
      
      #accmatequal[j, i] = Accuracy(all_true, equal_weights, test_labels, n.sp)
      #RSSequal[j, i] = RSS(pred.equal, test_labels)
      #meanRSSequal[j, i] = RSS(pred.equal, test_labels)/length(test_labels)
      
      accmatCF[j, i] = Accuracy(all_true, CF_weights, test_labels, n.sp)
      RSSCF[j, i] = RSS(pred.CF, test_labels)
      meanRSSCF[j, i] = RSS(pred.CF, test_labels)/length(test_labels)
      
      accmatkps[j, i] = Accuracy(all_true, kps_weights, test_labels, n.sp)
      RSSkps[j, i] = RSS(pred.kps, test_labels)
      meanRSSkps[j, i] = RSS(pred.kps, test_labels)/length(test_labels)
      
      #--
      if(is.null(cov.bias)){
        pred.knn = predict(simknn$fit.final, locations = sp_int_im)
      }else{
        pred.knn = predict(simknn$fit.final, covariates = pred.list, locations = sp_int_im)
      }
      
      sp.predlist.knn = list()
      for (l in 1:n.sp) {
        sp.predlist.knn[[l]] = as.vector(t(pred.knn[[l]]$v))
        
        IMSEknn[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.knn[[l]]/(simknn$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        sumcorknn1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.knn[[l]], method="pearson"))
        sumcorknn2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.knn[[l]], method="kendall"))
        
      }
      
      #--
      if(is.null(cov.bias)){
        pred.kmeans = predict(simkmeans$fit.final, locations = sp_int_im)
      }else{
        pred.kmeans = predict(simkmeans$fit.final, covariates = pred.list, locations = sp_int_im)
      }
      
      sp.predlist.kmeans = list()
      for (l in 1:n.sp) {
        sp.predlist.kmeans[[l]] = as.vector(t(pred.kmeans[[l]]$v))
        
        IMSEkmeans[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.kmeans[[l]]/(simkmeans$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        sumcorkmeans1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.kmeans[[l]], method="pearson"))
        sumcorkmeans2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.kmeans[[l]], method="kendall"))
        
      }
      
      #--
      if(is.null(cov.bias)){
        pred.random  = predict(simrandom$fit.final, locations = sp_int_im)
      }else{
        pred.random = predict(simrandom$fit.final, covariates = pred.list, locations = sp_int_im)
      }
      
      sp.predlist.rand = list()
      for (l in 1:n.sp) {
        sp.predlist.rand[[l]] = as.vector(t(pred.random[[l]]$v))
        
        IMSErand[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.rand[[l]]/(simrandom$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        sumcorrand1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.rand[[l]], method="pearson"))
        sumcorrand2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.rand[[l]], method="kendall"))
        
      }
      
      #--
      # if(is.null(cov.bias)){
      #   pred.equal = predict(simequal$fit.final, locations = sp1_int_im)
      # }else{
      #   pred.equal = predict(simequal$fit.final, covariates = pred.list, locations = sp1_int_im)
      # }
      # 
      # sp.predlist.equal = list()
      # for (l in 1:n.sp) {
      #   sp.predlist.equal[[l]] = as.vector(t(pred.equal[[l]]$v))
      #   
      #   IMSEequal[j, i] = sum(IMSE(sp_int.list[[l]], sp.predlist.equal[[l]]))
      #   sumcorequal1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.equal[[l]], method="pearson"))
      #   sumcorequal2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.equal[[l]], method="kendall"))
      #   
      # }
      
      #--
      if(is.null(cov.bias)){
        pred.CF  = predict(simCF$fit.final, locations = sp_int_im)
      }else{
        pred.CF = predict(simCF$fit.final, covariates = pred.list, 
                          locations = sp_int_im)
      }
      
      sp.predlist.CF= list()
      for (l in 1:n.sp) {
        sp.predlist.CF[[l]] = as.vector(t(pred.CF[[l]]$v))
        
        IMSECF[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.CF[[l]]/(simCF$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        sumcorCF1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.CF[[l]], 
                               method="pearson"))
        sumcorCF2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.CF[[l]], 
                               method="kendall"))
        
      }
      
      #--
      if(is.null(cov.bias)){
        pred.kps = predict(simkps$fit.final, locations = sp_int_im)
      }else{
        pred.kps = predict(simkps$fit.final, covariates = pred.list, 
                           locations = sp_int_im)
      }
      
      sp.predlist.kps = list()
      for (l in 1:n.sp) {
        sp.predlist.kps[[l]] = as.vector(t(pred.kps[[l]]$v))
        
        IMSEkps[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.kps[[l]]/(simkps$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        sumcorkps1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.kps[[l]], 
                                method="pearson"))
        sumcorkps2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.kps[[l]], 
                                method="kendall"))
        
      }
      
      # for intensity plots
      knnpred[,,j,i] = as.matrix(unlist(sp.predlist.knn))
      kmeanspred[,,j,i] =as.matrix(unlist(sp.predlist.kmeans))
      randpred[,,j,i] = as.matrix(unlist(sp.predlist.rand))
      #equalpred[,,j,i] = as.matrix(unlist(sp.predlist.equal))
      CFpred[,,j,i] = as.matrix(unlist(sp.predlist.CF))
      kpspred[,,j,i] = as.matrix(unlist(sp.predlist.kps))
      
      # for coefficients
      #coef.knn.mat[,,j,i] = as.matrix(simknn$fit.final$coef)
      #coef.kmeans.mat[,,j,i] = as.matrix(simkmeans$fit.final$coef)
      #coef.rand.mat[,,j,i] = as.matrix(simrandom$fit.final$coef)
      #coef.eq.mat[,,j,i] = as.matrix(simequal$fit.final$coef)
      #coef.CF.mat[,,j,i] = as.matrix(simCF$fit.final$coef)
      #coef.kps.mat[,,j,i] = as.matrix(simkps$fit.final$coef)      
      
      #---
      # ppmLoopEngine
      ###---
      
      simLoopA = ppmLoopEngine(datappp, all_test, n.sp, addpt = "LoopA", quads.= quads.win,
                               ppmform=ppmform, delta_max=NULL, delta_min=NULL, delta_step=NULL, win. = win, num.add = NULL,
                               cov.list.=cov.list, cov.bias=NULL, kVal =NULL, kAreaInt=NULL, maxit = 50,
                               tol=0.000001, verbose = TRUE, plots = FALSE)
      
      simLoopT = ppmLoopEngine(datappp, all_test, n.sp, addpt = "LoopT", quads.= quads.win,
                               ppmform= ppmform, delta_max=delta_max, delta_min=delta_min, delta_step=delta_step, win. = win, num.add = NULL,
                               cov.list.=cov.list, cov.bias=NULL, kVal =NULL, kAreaInt=NULL, maxit = 50,
                               tol=0.000001, verbose = TRUE, plots = FALSE)
      
      simLoopE = ppmLoopEngine(datappp, all_test, n.sp, addpt = "LoopE", quads.= quads.win,
                               ppmform= ppmform, delta_max=NULL, delta_min=NULL, delta_step=NULL, win. = win, num.add = num.add,
                               cov.list.=cov.list, cov.bias=NULL, kVal =NULL, kAreaInt=NULL, maxit = 50,
                               tol=0.000001, verbose = TRUE, plots = FALSE)
      
      
      # for performance measures
      
      LoopA_weights = simLoopA$New_weights
      pred.LoopA    = LoopA_weights
      colnames(pred.LoopA)= c(unique(test_labels))
      W.LA[[i]][[j]] = LoopA_weights
      
      LoopT_weights = simLoopT$New_weights
      pred.LoopT    = LoopT_weights
      colnames(pred.LoopT)= c(unique(test_labels))
      W.LT[[i]][[j]] = LoopT_weights
      
      LoopE_weights = simLoopE$New_weights
      pred.LoopE    = LoopE_weights
      colnames(pred.LoopE)= c(unique(test_labels))
      W.LE[[i]][[j]] = LoopE_weights
      
      #
      accmatLoopA[j, i] = Accuracy(all_true, LoopA_weights, test_labels, n.sp)
      RSSLoopA[j, i] = RSS(pred.LoopA, test_labels)
      meanRSSLoopA[j, i] = RSS(pred.LoopA, test_labels)/length(test_labels)
      
      accmatLoopT[j, i] = Accuracy(all_true, LoopT_weights, test_labels, n.sp)
      RSSLoopT[j, i] = RSS(pred.LoopT, test_labels)
      meanRSSLoopT[j, i] = RSS(pred.LoopT, test_labels)/length(test_labels)
      
      accmatLoopE[j, i] = Accuracy(all_true, LoopE_weights, test_labels, n.sp)
      RSSLoopE[j, i] = RSS(pred.LoopE, test_labels)
      meanRSSLoopE[j, i] = RSS(pred.LoopE, test_labels)/length(test_labels)
      
      #--
      pr_quad_ppmlist.LA = pr_quad_ppmlist.LT = pr_quad_ppmlist.LE = list()
      for (l in 1:n.sp) {
        if(is.null(cov.bias)){
          pr_quad_ppmlist.LA[[l]] = predict(simLoopA$ppm_list[[l]], locations = quads.win)
        }else{
          pr_quad_ppmlist.LA[[l]] = predict(simLoopA$ppm_list[[l]], covariates = pred.list, locations = quads.win)
        }
        
        if(is.null(cov.bias)){
          pr_quad_ppmlist.LT[[l]] = predict(simLoopT$ppm_list[[l]], locations = quads.win)
        }else{
          pr_quad_ppmlist.LT[[l]] = predict(simLoopT$ppm_list[[l]], covariates = pred.list, locations = quads.win)
        }
        
        if(is.null(cov.bias)){
          pr_quad_ppmlist.LE[[l]] = predict(simLoopE$ppm_list[[l]], locations = quads.win)
        }else{
          pr_quad_ppmlist.LE[[l]] = predict(simLoopE$ppm_list[[l]], covariates = pred.list, locations = quads.win)
        }
      }
      
      # for intensity plots
      LoopApred[,,j,i] = matrix(unlist(pr_quad_ppmlist.LA),
                                nrow=length(pr_quad_ppmlist.LA[[1]]), byrow=F)
      
      LoopTpred[,,j,i] = matrix(unlist(pr_quad_ppmlist.LT),
                                nrow=length(pr_quad_ppmlist.LT[[1]]), byrow=F)
      
      LoopEpred[,,j,i] = matrix(unlist(pr_quad_ppmlist.LE),
                                nrow=length(pr_quad_ppmlist.LE[[1]]), byrow=F)
      
      
      coef.LAvec = coef.LTvec = coef.LEvec = list()
      for (l in 1:n.sp) {
        # for coefficients
        coef.LAvec[[l]] = as.vector(unlist(simLoopA$ppm_list[[l]]$coef))
        coef.LTvec[[l]] = as.vector(unlist(simLoopT$ppm_list[[l]]$coef))
        coef.LEvec[[l]] = as.vector(unlist(simLoopE$ppm_list[[l]]$coef))
        
        l=l+1
      }
      
      coef.LA.mat[,,j,i] = as.vector(unlist(t(coef.LAvec)))
      coef.LT.mat[,,j,i] = as.vector(unlist(t(coef.LTvec)))
      coef.LE.mat[,,j,i] = as.vector(unlist(t(coef.LEvec)))
      
      sp.predlist.LA = sp.predlist.LT = sp.predlist.LE = list()
      for (l in 1:n.sp) {
        sp.predlist.LA[[l]] = as.vector(t(pr_quad_ppmlist.LA[[l]]))
        IMSELoopA[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.LA[[l]]/(length(simLoopA$sp_aug.list[[l]]$X)))*datappp$n/n.sp))
        sumcorLoopA1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.LA[[l]]))
        sumcorLoopA2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.LA[[l]]))
        
        sp.predlist.LT[[l]] = as.vector(t(pr_quad_ppmlist.LT[[l]]))
        IMSELoopT[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.LT[[l]]/(length(simLoopT$sp_aug.list[[l]]$X)))*datappp$n/n.sp))
        sumcorLoopT1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.LT[[l]]))
        sumcorLoopT2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.LT[[l]]))
        
        sp.predlist.LE[[l]] = as.vector(t(pr_quad_ppmlist.LE[[l]]))
        IMSELoopE[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predlist.LE[[l]]/(length(simLoopE$sp_aug.list[[l]]$X)))*datappp$n/n.sp))
        sumcorLoopE1[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.LE[[l]]))
        sumcorLoopE2[j, i] = sum(corint(sp_int.list[[l]], sp.predlist.LE[[l]]))
        
      }
      
      ####
      # Individual PPMs
      ###---
      
      datamarks = marks(datappp)
      uniquemarks = unique(datamarks)
      colmarks = uniquemarks
      unknown = datamarks == "Unknown"
      names.mark  = colmarks[-which(colmarks == "Unknown")]
      cov.list. = cov.list 
	  
      Qind = ppm_spind = pr_sp_ppmind = list()
      for (l in 1:n.sp) {
        #specie separetely
        Qind[[l]] = quadscheme(data = sp_sub.list[[l]], dummy = quads.win, method = "grid", ntile = c(101, 101), npix = c(101, 101))
        
        # Fit Poisson PPMs
        if(is.null(kAreaInt)){
          ppm_spind[[l]]  = ppm(Qind[[l]], trend = ppmform, covariates = cov.list.)
          
        }else{
          ppm_spind[[l]]  = ppm(Qind[[l]], trend = ppmform, covariates = cov.list., AreaInter(kAreaInt))
        }
        
        # Predict intensity at locations where labels are hidden
        ## former label from sp1 tested for the 3PPMs
        if(is.null(cov.bias)){
          pr_sp_ppmind[[l]] = predict(ppm_spind[[l]], locations = all_test)
        }else{
          pr_sp_ppmind[[l]] = predict(ppm_spind[[l]], covariates = pred.list, locations = all_test)
        }
      }
      
      
      sp_predsind = data.frame(matrix(unlist(pr_sp_ppmind),
                                      nrow=length(pr_sp_ppmind[[1]]), byrow=F))
      ind_wts  = sp_predsind/apply(sp_predsind, 1, sum)
      max_predind = apply(sp_predsind, 1, which.max)
      max_predind.vec = as.vector(table(max_predind))
      
      # calculate accuracy in the same way that we did for the mixture models
      check.id = as.vector(max_predind)
      check.data = as.data.frame(cbind(test_labels, check.id))
      
      # to deal with some labels not be present in some iterations
      for (l in 1:n.sp) {
        if (anyNA(sum(check.data$test_labels == paste("Sp", l, sep="")))){
          allvec = as.vector(seq(from=1, to=n.sp, by=1))
          levels(check.data$test_labels) = allvec[which(allvec!=l)]
        }else{
          allvec = as.vector(seq(from=1, to=n.sp, by=1))
          levels(check.data$test_labels) = allvec
        }
      }
      
      testlab = as.data.frame(check.data$test_labels)
      checklab = as.data.frame(check.data$check.id)
      
      checkD = as.data.frame(check.data)
      
      m.acc.indiv = createConfusionMatrix(checkD$test_labels, checkD$check.id)
      
      n.check = sum(m.acc.indiv) # number of instances
      diag.check = diag(m.acc.indiv) # number of correctly classified instances per class
      
      # accuracy measure
      accmatindiv[j, i] = sum(diag.check) / n.check
      
      
      # RSS measure
      # Assign weights
      all_predsind = sp_predsind
      testind_wts  = all_predsind/apply(all_predsind, 1, sum)
      colnames(testind_wts) = names.mark
      W.ind[[i]][[j]] = testind_wts
      
      RSSindiv[j, i] = RSS(testind_wts, test_labels)
      meanRSSindiv[j, i] = RSS(testind_wts, test_labels)/length(test_labels)
      
      pr_quad_ppmindlist = list()
      for (l in 1:n.sp) {
        if(is.null(cov.bias)){
          pr_quad_ppmindlist[[l]] = predict(ppm_spind[[l]], locations = quads.win)
        }else{
          pr_quad_ppmindlist[[l]] = predict(ppm_spind[[l]], covariates = pred.list, locations = quads.win)
        }
        
      }
      
      # for intensity plots
      indivpred[,,j,i] = matrix(unlist(pr_quad_ppmindlist),
                                nrow=length(pr_quad_ppmindlist[[1]]), byrow=F)
      
      coef.indvec = list()
      for (l in 1:n.sp) {
        # for coefficients
        coef.indvec[[l]] = as.vector(unlist(ppm_spind[[l]]$coef))
        
        l=l+1
      }
      
      coef.ind.mat[,,j,i] = as.vector(unlist(t(coef.indvec)))
      
      sp.predindlist = list()
      for (l in 1:n.sp) {
        sp.predindlist[[l]] = as.vector(t(pr_quad_ppmindlist[[l]]))
        IMSEindiv[j, i] = sum(IMSE(sp_int.list[[l]], (sp.predindlist[[l]]/(length(sp_sub.list[[l]]$x)))*datappp$n/n.sp))        
        sumcorindiv1[j, i] = sum(corint(sp_int.list[[l]], sp.predindlist[[l]]))
        sumcorindiv2[j, i] = sum(corint(sp_int.list[[l]], sp.predindlist[[l]]))
        
      }
      
      # Calculate Standard errors
      if(length(hidepct)==1){
       #at locations 
        un.knn1 <- predict(simknn$fit.final[[1]], locations=simknn$sp_aug.list[[1]], se=TRUE)
        un.kmeans1 <- predict(simkmeans$fit.final[[1]], locations=simkmeans$sp_aug.list[[1]], se=TRUE)
        un.rand1 <- predict(simrandom$fit.final[[1]], locations=simrandom$sp_aug.list[[1]], se=TRUE)
        un.CF1 <- predict(simCF$fit.final[[1]], locations=simCF$sp_aug.list[[1]], se=TRUE)
        un.kps1 <- predict(simkps$fit.final[[1]], locations=simkps$sp_aug.list[[1]], se=TRUE)
        un.LA1 <- predict(simLoopA$ppm_list[[1]], locations=simLoopA$sp_aug_ppp.list[[1]], se=TRUE)
        un.LT1 <- predict(simLoopT$ppm_list[[1]], locations=simLoopT$sp_aug_ppp.list[[1]], se=TRUE)
        un.LE1 <- predict(simLoopE$ppm_list[[1]], locations=simLoopE$sp_aug_ppp.list[[1]], se=TRUE)
        un.indiv1 <- predict(ppm_spind[[1]], locations=datappp[datappp$marks=="Sp1"], se=TRUE)
        
        se.knn1[[j]] <- un.knn1$se/simknn$sp_aug.list[[1]]$n
        se.kmeans1[[j]] <- un.kmeans1$se/simkmeans$sp_aug.list[[1]]$n
        se.rand1[[j]] <- un.rand1$se/simrandom$sp_aug.list[[1]]$n
        se.CF1[[j]] <- un.CF1$se/simCF$sp_aug.list[[1]]$n
        se.kps1[[j]] <- un.kps1$se/simkps$sp_aug.list[[1]]$n
        se.LA1[[j]] <- un.LA1$se/simLoopA$sp_aug_ppp.list[[1]]$n
        se.LT1[[j]] <- un.LT1$se/simLoopT$sp_aug_ppp.list[[1]]$n
        se.LE1[[j]] <- un.LE1$se/simLoopE$sp_aug_ppp.list[[1]]$n
        se.indiv1[[j]] <- un.indiv1$se/datappp[datappp$marks=="Sp1"]$n
        
        #all
        w.knn1 <- predict(simknn$fit.final[[1]], se=TRUE)
        w.kmeans1 <- predict(simkmeans$fit.final[[1]], se=TRUE)
        w.rand1 <- predict(simrandom$fit.final[[1]], se=TRUE)
        w.CF1 <- predict(simCF$fit.final[[1]], se=TRUE)
        w.kps1 <- predict(simkps$fit.final[[1]], se=TRUE)
        w.LA1 <- predict(simLoopA$ppm_list[[1]], se=TRUE)
        w.LT1 <- predict(simLoopT$ppm_list[[1]], se=TRUE)
        w.LE1 <- predict(simLoopE$ppm_list[[1]], se=TRUE)
        w.indiv1 <- predict(ppm_spind[[1]], se=TRUE)
        
        sew.knn1[[j]] <- w.knn1$se
        sew.kmeans1[[j]] <- w.kmeans1$se
        sew.rand1[[j]] <- w.rand1$se
        sew.CF1[[j]] <- w.CF1$se
        sew.kps1[[j]] <- w.kps1$se
        sew.LA1[[j]] <- w.LA1$se
        sew.LT1[[j]] <- w.LT1$se
        sew.LE1[[j]] <- w.LE1$se
        sew.indiv1[[j]] <- w.indiv1$se
      }else{
        un.knn1 <- predict(simknn$fit.final[[1]], locations=simknn$sp_aug.list[[1]], se=TRUE)
        un.kmeans1 <- predict(simkmeans$fit.final[[1]], locations=simkmeans$sp_aug.list[[1]], se=TRUE)
        un.rand1 <- predict(simrandom$fit.final[[1]], locations=simrandom$sp_aug.list[[1]], se=TRUE)
        un.CF1 <- predict(simCF$fit.final[[1]], locations=simCF$sp_aug.list[[1]], se=TRUE)
        un.kps1 <- predict(simkps$fit.final[[1]], locations=simkps$sp_aug.list[[1]], se=TRUE)
        un.LA1 <- predict(simLoopA$ppm_list[[1]], locations=simLoopA$sp_aug_ppp.list[[1]], se=TRUE)
        un.LT1 <- predict(simLoopT$ppm_list[[1]], locations=simLoopT$sp_aug_ppp.list[[1]], se=TRUE)
        un.LE1 <- predict(simLoopE$ppm_list[[1]], locations=simLoopE$sp_aug_ppp.list[[1]], se=TRUE)
        un.indiv1 <- predict(ppm_spind[[1]], locations=datappp[datappp$marks=="Sp1"], se=TRUE)
        
        se.knn1[[i]][[j]] <- un.knn1$se/simknn$sp_aug.list[[1]]$n
        se.kmeans1[[i]][[j]] <- un.kmeans1$se/simkmeans$sp_aug.list[[1]]$n
        se.rand1[[i]][[j]] <- un.rand1$se/simrandom$sp_aug.list[[1]]$n
        se.CF1[[i]][[j]] <- un.CF1$se/simCF$sp_aug.list[[1]]$n
        se.kps1[[i]][[j]] <- un.kps1$se/simkps$sp_aug.list[[1]]$n
        se.LA1[[i]][[j]] <- un.LA1$se/simLoopA$sp_aug_ppp.list[[1]]$n
        se.LT1[[i]][[j]] <- un.LT1$se/simLoopT$sp_aug_ppp.list[[1]]$n
        se.LE1[[i]][[j]] <- un.LE1$se/simLoopE$sp_aug_ppp.list[[1]]$n
        se.indiv1[[i]][[j]] <- un.indiv1$se/datappp[datappp$marks=="Sp1"]$n
        
        #all
        w.knn1 <- predict(simknn$fit.final[[1]], se=TRUE)
        w.kmeans1 <- predict(simkmeans$fit.final[[1]], se=TRUE)
        w.rand1 <- predict(simrandom$fit.final[[1]], se=TRUE)
        w.CF1 <- predict(simCF$fit.final[[1]], se=TRUE)
        w.kps1 <- predict(simkps$fit.final[[1]], se=TRUE)
        w.LA1 <- predict(simLoopA$ppm_list[[1]], se=TRUE)
        w.LT1 <- predict(simLoopT$ppm_list[[1]], se=TRUE)
        w.LE1 <- predict(simLoopE$ppm_list[[1]], se=TRUE)
        w.indiv1 <- predict(ppm_spind[[1]], se=TRUE)
        
        sew.knn1[[i]][[j]] <- w.knn1$se
        sew.kmeans1[[i]][[j]] <- w.kmeans1$se
        sew.rand1[[i]][[j]] <- w.rand1$se
        sew.CF1[[i]][[j]] <- w.CF1$se
        sew.kps1[[i]][[j]] <- w.kps1$se
        sew.LA1[[i]][[j]] <- w.LA1$se
        sew.LT1[[i]][[j]] <- w.LT1$se
        sew.LE1[[i]][[j]] <- w.LE1$se
        sew.indiv1[[i]][[j]] <- w.indiv1$se
      }
      
      if(length(hidepct)==1){
       #at locations 
        un.knn2 <- predict(simknn$fit.final[[2]], locations=simknn$sp_aug.list[[2]], se=TRUE)
        un.kmeans2 <- predict(simkmeans$fit.final[[2]], locations=simkmeans$sp_aug.list[[2]], se=TRUE)
        un.rand2 <- predict(simrandom$fit.final[[2]], locations=simrandom$sp_aug.list[[2]], se=TRUE)
        un.CF2 <- predict(simCF$fit.final[[2]], locations=simCF$sp_aug.list[[2]], se=TRUE)
        un.kps2 <- predict(simkps$fit.final[[2]], locations=simkps$sp_aug.list[[2]], se=TRUE)
        un.LA2 <- predict(simLoopA$ppm_list[[2]], locations=simLoopA$sp_aug_ppp.list[[2]], se=TRUE)
        un.LT2 <- predict(simLoopT$ppm_list[[2]], locations=simLoopT$sp_aug_ppp.list[[2]], se=TRUE)
        un.LE2 <- predict(simLoopE$ppm_list[[2]], locations=simLoopE$sp_aug_ppp.list[[2]], se=TRUE)
        un.indiv2 <- predict(ppm_spind[[2]], locations=datappp[datappp$marks=="Sp2"], se=TRUE)
        
        se.knn2[[j]] <- un.knn2$se/simknn$sp_aug.list[[2]]$n
        se.kmeans2[[j]] <- un.kmeans2$se/simkmeans$sp_aug.list[[2]]$n
        se.rand2[[j]] <- un.rand2$se/simrandom$sp_aug.list[[2]]$n
        se.CF2[[j]] <- un.CF2$se/simCF$sp_aug.list[[2]]$n
        se.kps2[[j]] <- un.kps2$se/simkps$sp_aug.list[[2]]$n
        se.LA2[[j]] <- un.LA2$se/simLoopA$sp_aug_ppp.list[[2]]$n
        se.LT2[[j]] <- un.LT2$se/simLoopT$sp_aug_ppp.list[[2]]$n
        se.LE2[[j]] <- un.LE2$se/simLoopE$sp_aug_ppp.list[[2]]$n
        se.indiv2[[j]] <- un.indiv2$se/datappp[datappp$marks=="Sp2"]$n
        
        #all
        w.knn2 <- predict(simknn$fit.final[[2]], se=TRUE)
        w.kmeans2 <- predict(simkmeans$fit.final[[2]], se=TRUE)
        w.rand2 <- predict(simrandom$fit.final[[2]], se=TRUE)
        w.CF2 <- predict(simCF$fit.final[[2]], se=TRUE)
        w.kps2 <- predict(simkps$fit.final[[2]], se=TRUE)
        w.LA2 <- predict(simLoopA$ppm_list[[2]], se=TRUE)
        w.LT2 <- predict(simLoopT$ppm_list[[2]], se=TRUE)
        w.LE2 <- predict(simLoopE$ppm_list[[2]], se=TRUE)
        w.indiv2 <- predict(ppm_spind[[2]], se=TRUE)
        
        sew.knn2[[j]] <- w.knn2$se
        sew.kmeans2[[j]] <- w.kmeans2$se
        sew.rand2[[j]] <- w.rand2$se
        sew.CF2[[j]] <- w.CF2$se
        sew.kps2[[j]] <- w.kps2$se
        sew.LA2[[j]] <- w.LA2$se
        sew.LT2[[j]] <- w.LT2$se
        sew.LE2[[j]] <- w.LE2$se
        sew.indiv2[[j]] <- w.indiv2$se
      }else{
        un.knn2 <- predict(simknn$fit.final[[2]], locations=simknn$sp_aug.list[[2]], se=TRUE)
        un.kmeans2 <- predict(simkmeans$fit.final[[2]], locations=simkmeans$sp_aug.list[[2]], se=TRUE)
        un.rand2 <- predict(simrandom$fit.final[[2]], locations=simrandom$sp_aug.list[[2]], se=TRUE)
        un.CF2 <- predict(simCF$fit.final[[2]], locations=simCF$sp_aug.list[[2]], se=TRUE)
        un.kps2 <- predict(simkps$fit.final[[2]], locations=simkps$sp_aug.list[[2]], se=TRUE)
        un.LA2 <- predict(simLoopA$ppm_list[[2]], locations=simLoopA$sp_aug_ppp.list[[2]], se=TRUE)
        un.LT2 <- predict(simLoopT$ppm_list[[2]], locations=simLoopT$sp_aug_ppp.list[[2]], se=TRUE)
        un.LE2 <- predict(simLoopE$ppm_list[[2]], locations=simLoopE$sp_aug_ppp.list[[2]], se=TRUE)
        un.indiv2 <- predict(ppm_spind[[2]], locations=datappp[datappp$marks=="Sp2"], se=TRUE)
        
        se.knn2[[i]][[j]] <- un.knn2$se/simknn$sp_aug.list[[2]]$n
        se.kmeans2[[i]][[j]] <- un.kmeans2$se/simkmeans$sp_aug.list[[2]]$n
        se.rand2[[i]][[j]] <- un.rand2$se/simrandom$sp_aug.list[[2]]$n
        se.CF2[[i]][[j]] <- un.CF2$se/simCF$sp_aug.list[[2]]$n
        se.kps2[[i]][[j]] <- un.kps2$se/simkps$sp_aug.list[[2]]$n
        se.LA2[[i]][[j]] <- un.LA2$se/simLoopA$sp_aug_ppp.list[[2]]$n
        se.LT2[[i]][[j]] <- un.LT2$se/simLoopT$sp_aug_ppp.list[[2]]$n
        se.LE2[[i]][[j]] <- un.LE2$se/simLoopE$sp_aug_ppp.list[[2]]$n
        se.indiv2[[i]][[j]] <- un.indiv2$se/datappp[datappp$marks=="Sp2"]$n
        
        #all
        w.knn2 <- predict(simknn$fit.final[[2]], se=TRUE)
        w.kmeans2 <- predict(simkmeans$fit.final[[2]], se=TRUE)
        w.rand2 <- predict(simrandom$fit.final[[2]], se=TRUE)
        w.CF2 <- predict(simCF$fit.final[[2]], se=TRUE)
        w.kps2 <- predict(simkps$fit.final[[2]], se=TRUE)
        w.LA2 <- predict(simLoopA$ppm_list[[2]], se=TRUE)
        w.LT2 <- predict(simLoopT$ppm_list[[2]], se=TRUE)
        w.LE2 <- predict(simLoopE$ppm_list[[2]], se=TRUE)
        w.indiv2 <- predict(ppm_spind[[2]], se=TRUE)
        
        sew.knn2[[i]][[j]] <- w.knn2$se
        sew.kmeans2[[i]][[j]] <- w.kmeans2$se
        sew.rand2[[i]][[j]] <- w.rand2$se
        sew.CF2[[i]][[j]] <- w.CF2$se
        sew.kps2[[i]][[j]] <- w.kps2$se
        sew.LA2[[i]][[j]] <- w.LA2$se
        sew.LT2[[i]][[j]] <- w.LT2$se
        sew.LE2[[i]][[j]] <- w.LE2$se
        sew.indiv2[[i]][[j]] <- w.indiv2$se
      }
      
	  if(length(hidepct)==2){
       #at locations 
        un.knn3 <- predict(simknn$fit.final[[3]], locations=simknn$sp_aug.list[[3]], se=TRUE)
        un.kmeans3 <- predict(simkmeans$fit.final[[3]], locations=simkmeans$sp_aug.list[[3]], se=TRUE)
        un.rand3 <- predict(simrandom$fit.final[[3]], locations=simrandom$sp_aug.list[[3]], se=TRUE)
        un.CF3 <- predict(simCF$fit.final[[3]], locations=simCF$sp_aug.list[[3]], se=TRUE)
        un.kps3 <- predict(simkps$fit.final[[3]], locations=simkps$sp_aug.list[[3]], se=TRUE)
        un.LA3 <- predict(simLoopA$ppm_list[[3]], locations=simLoopA$sp_aug_ppp.list[[3]], se=TRUE)
        un.LT3 <- predict(simLoopT$ppm_list[[3]], locations=simLoopT$sp_aug_ppp.list[[3]], se=TRUE)
        un.LE3 <- predict(simLoopE$ppm_list[[3]], locations=simLoopE$sp_aug_ppp.list[[3]], se=TRUE)
        un.indiv3 <- predict(ppm_spind[[3]], locations=datappp[datappp$marks=="Sp3"], se=TRUE)
        
        se.knn3[[j]] <- un.knn3$se/simknn$sp_aug.list[[3]]$n
        se.kmeans3[[j]] <- un.kmeans3$se/simkmeans$sp_aug.list[[3]]$n
        se.rand3[[j]] <- un.rand3$se/simrandom$sp_aug.list[[3]]$n
        se.CF3[[j]] <- un.CF3$se/simCF$sp_aug.list[[3]]$n
        se.kps3[[j]] <- un.kps3$se/simkps$sp_aug.list[[3]]$n
        se.LA3[[j]] <- un.LA3$se/simLoopA$sp_aug_ppp.list[[3]]$n
        se.LT3[[j]] <- un.LT3$se/simLoopT$sp_aug_ppp.list[[3]]$n
        se.LE3[[j]] <- un.LE3$se/simLoopE$sp_aug_ppp.list[[3]]$n
        se.indiv3[[j]] <- un.indiv3$se/datappp[datappp$marks=="Sp3"]$n
        
        #all
        w.knn3 <- predict(simknn$fit.final[[3]], se=TRUE)
        w.kmeans3 <- predict(simkmeans$fit.final[[3]], se=TRUE)
        w.rand3 <- predict(simrandom$fit.final[[3]], se=TRUE)
        w.CF3 <- predict(simCF$fit.final[[3]], se=TRUE)
        w.kps3 <- predict(simkps$fit.final[[3]], se=TRUE)
        w.LA3 <- predict(simLoopA$ppm_list[[3]], se=TRUE)
        w.LT3 <- predict(simLoopT$ppm_list[[3]], se=TRUE)
        w.LE3 <- predict(simLoopE$ppm_list[[3]], se=TRUE)
        w.indiv3 <- predict(ppm_spind[[3]], se=TRUE)
        
        sew.knn3[[j]] <- w.knn3$se
        sew.kmeans3[[j]] <- w.kmeans3$se
        sew.rand3[[j]] <- w.rand3$se
        sew.CF3[[j]] <- w.CF3$se
        sew.kps3[[j]] <- w.kps3$se
        sew.LA3[[j]] <- w.LA3$se
        sew.LT3[[j]] <- w.LT3$se
        sew.LE3[[j]] <- w.LE3$se
        sew.indiv3[[j]] <- w.indiv3$se
      }else{
        un.knn3 <- predict(simknn$fit.final[[3]], locations=simknn$sp_aug.list[[3]], se=TRUE)
        un.kmeans3 <- predict(simkmeans$fit.final[[3]], locations=simkmeans$sp_aug.list[[3]], se=TRUE)
        un.rand3 <- predict(simrandom$fit.final[[3]], locations=simrandom$sp_aug.list[[3]], se=TRUE)
        un.CF3 <- predict(simCF$fit.final[[3]], locations=simCF$sp_aug.list[[3]], se=TRUE)
        un.kps3 <- predict(simkps$fit.final[[3]], locations=simkps$sp_aug.list[[3]], se=TRUE)
        un.LA3 <- predict(simLoopA$ppm_list[[3]], locations=simLoopA$sp_aug_ppp.list[[3]], se=TRUE)
        un.LT3 <- predict(simLoopT$ppm_list[[3]], locations=simLoopT$sp_aug_ppp.list[[3]], se=TRUE)
        un.LE3 <- predict(simLoopE$ppm_list[[3]], locations=simLoopE$sp_aug_ppp.list[[3]], se=TRUE)
        un.indiv3 <- predict(ppm_spind[[3]], locations=datappp[datappp$marks=="Sp3"], se=TRUE)
        
        se.knn3[[i]][[j]] <- un.knn3$se/simknn$sp_aug.list[[3]]$n
        se.kmeans3[[i]][[j]] <- un.kmeans3$se/simkmeans$sp_aug.list[[3]]$n
        se.rand3[[i]][[j]] <- un.rand3$se/simrandom$sp_aug.list[[3]]$n
        se.CF3[[i]][[j]] <- un.CF3$se/simCF$sp_aug.list[[3]]$n
        se.kps3[[i]][[j]] <- un.kps3$se/simkps$sp_aug.list[[3]]$n
        se.LA3[[i]][[j]] <- un.LA3$se/simLoopA$sp_aug_ppp.list[[3]]$n
        se.LT3[[i]][[j]] <- un.LT3$se/simLoopT$sp_aug_ppp.list[[3]]$n
        se.LE3[[i]][[j]] <- un.LE3$se/simLoopE$sp_aug_ppp.list[[3]]$n
        se.indiv3[[i]][[j]] <- un.indiv3$se/datappp[datappp$marks=="Sp3"]$n
		
        #all
        w.knn3 <- predict(simknn$fit.final[[3]], se=TRUE)
        w.kmeans3 <- predict(simkmeans$fit.final[[3]], se=TRUE)
        w.rand3 <- predict(simrandom$fit.final[[3]], se=TRUE)
        w.CF3 <- predict(simCF$fit.final[[3]], se=TRUE)
        w.kps3 <- predict(simkps$fit.final[[3]], se=TRUE)
        w.LA3 <- predict(simLoopA$ppm_list[[3]], se=TRUE)
        w.LT3 <- predict(simLoopT$ppm_list[[3]], se=TRUE)
        w.LE3 <- predict(simLoopE$ppm_list[[3]], se=TRUE)
        w.indiv3 <- predict(ppm_spind[[3]], se=TRUE)
        
        sew.knn3[[i]][[j]] <- w.knn3$se
        sew.kmeans3[[i]][[j]] <- w.kmeans3$se
        sew.rand3[[i]][[j]] <- w.rand3$se
        sew.CF3[[i]][[j]] <- w.CF3$se
        sew.kps3[[i]][[j]] <- w.kps3$se
        sew.LA3[[i]][[j]] <- w.LA3$se
        sew.LT3[[i]][[j]] <- w.LT3$se
        sew.LE3[[i]][[j]] <- w.LE3$se
        sew.indiv3[[i]][[j]] <- w.indiv3$se
      }
      cat(paste(i, j, "\n"))
      flush.console()
    }
  }
  return(list(RSSknn = RSSknn, meanRSSknn = meanRSSknn, IMSEknn = IMSEknn, sumcorknn2 = sumcorknn2, sumcorknn2 = sumcorknn2,
              RSSkmeans = RSSkmeans, meanRSSkmeans = meanRSSkmeans, IMSEkmeans = IMSEkmeans, sumcorkmeans2 = sumcorkmeans2, sumcorkmeans2 = sumcorkmeans2, 
              RSSrand = RSSrand, meanRSSrand = meanRSSrand, IMSErand = IMSErand, sumcorrand2 = sumcorrand2, sumcorrand2 = sumcorrand2,
              #RSSequal = RSSequal, meanRSSequal = meanRSSequal, IMSEequal = IMSEequal, sumcorequal2 = sumcorequal2, sumcorequal2 = sumcorequal2, 
              RSSkps = RSSkps, meanRSSkps = meanRSSkps, IMSEkps = IMSEkps, sumcorkps2 = sumcorkps2, sumcorkps2 = sumcorkps2,
              RSSCF = RSSCF, meanRSSCF = meanRSSCF, IMSECF = IMSECF, sumcorCF2 = sumcorCF2, sumcorCF2 = sumcorCF2,
              RSSindiv = RSSindiv, meanRSSindiv = meanRSSindiv, IMSEindiv = IMSEindiv, sumcorindiv2 = sumcorindiv2, sumcorindiv2 = sumcorindiv2, 
              RSSLoopA = RSSLoopA, meanRSSLoopA = meanRSSLoopA, IMSELoopA = IMSELoopA, sumcorLoopA2 = sumcorLoopA2, sumcorLoopA2 = sumcorLoopA2,
              RSSLoopT = RSSLoopT, meanRSSLoopT = meanRSSLoopT, IMSELoopT = IMSELoopT, sumcorLoopT2 = sumcorLoopT2, sumcorLoopT2 = sumcorLoopT2, 
              RSSLoopE = RSSLoopE, meanRSSLoopE = meanRSSLoopE, IMSELoopE = IMSELoopE, sumcorLoopE2 = sumcorLoopE2, sumcorLoopE2 = sumcorLoopE2, 
              accmatknn = accmatknn, accmatkmeans = accmatkmeans, accmatrand = accmatrand, #accmatequal = accmatequal, 
              accmatindiv = accmatindiv, accmatLoopA = accmatLoopA, accmatLoopT = accmatLoopT,accmatLoopE = accmatLoopE, accmatkps = accmatkps, accmatCF = accmatCF,
              knnpred = knnpred, kmeanspred = kmeanspred, randpred = randpred, #equalpred = equalpred, 
              LoopApred = LoopApred, LoopTpred = LoopTpred, LoopEpred = LoopEpred, kpspred = kpspred, CFpred = CFpred,
              indivpred = indivpred, coef.knn.mat = coef.knn.mat, coef.kmeans.mat = coef.kmeans.mat, 
              coef.rand.mat = coef.rand.mat, oef.LA.mat = coef.LA.mat,  
              coef.LT.mat = coef.LT.mat,  coef.LE.mat = coef.LE.mat,  coef.ind.mat = coef.ind.mat,  coef.kps.mat = coef.kps.mat,  coef.CF.mat = coef.CF.mat,#coef.eq.mat = coef.eq.mat, 
              W.knn = W.knn, W.kmeans = W.kmeans, W.rand = W.rand, W.kps = W.kps, W.CF = W.CF, W.ind = W.ind, W.LA = W.LA, W.LT = W.LT, W.LE = W.LE, sp.true = sp.true, se.knn2 = se.knn2, se.kmeans2 = se.kmeans2, se.rand2 = se.rand2, 
              se.CF1 = se.CF1, se.kps1 = se.kps1, se.LA1 = se.LA1, se.LT1 = se.LT1, se.LE1 = se.LE1, se.indiv1 = se.indiv1, se.knn1 = se.knn1, se.kmeans1 = se.kmeans1, se.rand1 = se.rand1, 
			  se.CF2 = se.CF2, se.kps2 = se.kps2, se.LA2 = se.LA2, se.LT2 = se.LT2, se.LE2 = se.LE2, se.indiv2 = se.indiv2, se.knn2 = se.knn2, se.kmeans2 = se.kmeans2, se.rand2 = se.rand2, 
              se.CF2 = se.CF2, se.kps2 = se.kps2, se.LA2 = se.LA2, se.LT2 = se.LT2, se.LE2 = se.LE2, se.indiv2 = se.indiv2, se.knn3 = se.knn3, se.kmeans3 = se.kmeans3, se.rand3 = se.rand3, 
              se.CF3 = se.CF3, se.kps3 = se.kps3, se.LA3 = se.LA3, se.LT3 = se.LT3, se.LE3 = se.LE3, se.indiv3 = se.indiv3, sew.knn1 = sew.knn1, sew.kmeans1 = sew.kmeans1, sew.rand1 = sew.rand1, 
              sew.CF1 = sew.CF1, sew.kps1 = sew.kps1, sew.LA1 = sew.LA1, sew.LT1 = sew.LT1, sew.LE1 = sew.LE1, sew.indiv1 = sew.indiv1, sew.knn2 = sew.knn2, sew.kmeans2 = sew.kmeans2, sew.rand2 = sew.rand2, 
			  sew.CF2 = sew.CF2, sew.kps2 = sew.kps2, sew.LA2 = sew.LA2, sew.LT2 = sew.LT2, sew.LE2 = sew.LE2, sew.indiv2 = sew.indiv2, sew.knn2 = sew.knn2, sew.kmeans2 = sew.kmeans2, sew.rand2 = sew.rand2, 
              sew.CF2 = sew.CF2, sew.kps2 = sew.kps2, sew.LA2 = sew.LA2, sew.LT2 = sew.LT2, sew.LE2 = sew.LE2, sew.indiv2 = sew.indiv2, sew.knn3 = sew.knn3, sew.kmeans3 = sew.kmeans3, sew.rand3 = sew.rand3, 
              sew.CF3 = sew.CF3, sew.kps3 = sew.kps3, sew.LA3 = sew.LA3, sew.LT3 = sew.LT3, sew.LE3 = sew.LE3, sew.indiv3 = sew.indiv3
              
  ))
}


QuickTest = Testsims(hidepct=c(0.2, 0.5, 0.8), n.sims=5, sp_sim.list, n.sp=3,
                     k = 1, cov.list=cov.list, cov.bias=NULL, kVal=NULL, kAreaInt=NULL,
					 delta_max=0.5, delta_min=0.1, delta_step =0.1, num.add = 5)