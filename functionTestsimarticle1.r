###########################################################################################
#
#                 Functions to use in the testsims function
#
###########################################################################################

# call the different function for the simulations


### Function ppmMixEngine
#---------------------------------------------------------------------

makeMask = function(Qppp)
{
  
  q_df = data.frame(x = Qppp$x, y = Qppp$y)
  ux = sort(unique(q_df$x))
  uy = sort(unique(q_df$y))
  nx = length(ux)
  ny = length(uy)
  
  col.ref = match(q_df$x, ux)
  row.ref = match(q_df$y, uy)
  
  all.vec          = rep(0, max(row.ref)*max(col.ref))
  vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
  all.vec[vec.ref] = 1
  mask.out         = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
  mask.out
}

#------------------------------------------------------------------------------------
#  								function ppmMixEngine
#------------------------------------------------------------------------------------

ppmMixEngine = function(datappp, quads., ppmform, all_true, all_test, n.sp, sp_int_im,
                        initweights = c("knn", "kps", "kmeans", "random", "CoinF"), 
                        k=k, ks=ks, nstart=nstart, cov.list., cov.bias=NULL, kVal = NULL, kAreaInt=NULL,
                        verbose = TRUE, tol = 0.001, maxit = 50, plots = FALSE)
{
  
  # define some objects
  datamarks = marks(datappp)
  uniquemarks = unique(datamarks)
  unknown = datamarks == "Unknown"
  names.mark = uniquemarks[uniquemarks != "Unknown"]
  nclust  = length(unique(datamarks)) - 1
  Qmask = makeMask(quads.)
  
  # separate the multi type point pattern
  splitppps = split(datappp, as.factor(marks(datappp)))
  for (i in 1:(nclust + 1))
  {
    assign(paste("ppp_", names(splitppps)[i], sep = ""),
           splitppps[[i]])
  }
  
  #1# Initialization of membership probabilities
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Set up otpion for initial weights  
  initweights <- match.arg(initweights)
  
  # for knn method
  if (initweights == "knn"){
    #1. Compute distances of unknown points to all known points
    if(k > all_true$n)
    {
      print("Impossible to consider so many points, the highest possible number will be used instead")
      k = all_true$n
    }
    
    knownmarks = datamarks[datamarks != "Unknown"] # marks of known points
    knownsp = sort(unique(knownmarks)) # names of known species
    n.sp = length(knownsp) # number of known species
    
    all_dists = as.matrix(dist(data.frame(datappp$x, datappp$y))) # all distances
    unknown_dists = all_dists[which(datamarks != "Unknown"), which(datamarks == "Unknown")] #distances between unkonwn points and known points
    
    #2. Compute unscaled membership weights for unknown species
    weight_first = matrix(NA, sum(datamarks == "Unknown"), n.sp) # initial unscaled weight
    
    for (i in 1:sum(datamarks == "Unknown"))
    {
      dist_i = unknown_dists[,i] # distances for ith unknown point
      shortest_k = sort(dist_i)[k] # find kth shortest distance
      for (j in 1:n.sp)
      {
        dist_i_sp_j = dist_i[knownmarks == knownsp[j]] # distances to species j
        within_shortest_k = dist_i_sp_j <= shortest_k # which are within the smallest k distances?
        shortdist_weights = 1/dist_i_sp_j*within_shortest_k # weight of 1/k for points within smallest k distances
        weight_first[i, j] = sum(shortdist_weights)
      }
    }
    
    #3. Rescale membership weights
    init_weight_unknown = prop.table(weight_first, 1) # rescale weights
    
    #4. Set up membership weights for points with known labels
    init_weight_known = matrix(0, length(knownmarks), n.sp) # weights for known species labels
    for (j in 1:n.sp)
    {
      init_weight_known[knownmarks == knownsp[j], j] = 1
    }
    
    #5. Combine
    init.weight = rbind(init_weight_known, init_weight_unknown) # combine
    colnames(init.weight) = knownsp
    
    iterweights = init.weight
    itermarks = datamarks
    itermarks = colnames(init.weight)[apply(init.weight, 1, which.max)]
    
  }
  
  # for kps method
  if (initweights == "kps"){
    if(max(ks) > min(table(all_true$marks)))
    {
      print("Impossible to consider so many points, the highest possible number will be used instead")
      ks = 1:min(table(all_true$marks))
    }
    
	# calculate the nearest neighbours distances
    nndists  = nndist(datappp, k=ks, by = as.factor(marks(datappp)))
    if(length(ks)>1){
      nndists.new = matrix(ncol=n.sp,nrow=datappp$n)
      nndists  = nndists[1:(n.sp*length(ks))]  # only consider known sp yet
      
      splitdf <- function(df, n) {
        indx <- matrix(seq_len(ncol(df)), ncol = n)
        lapply(seq_len(n), function(x) df[, indx[, x]])
      }
      
      nndist.sub = splitdf(nndists, n.sp) # split the dataset in different dataset according to the number of species
      
      for (s in 1:n.sp) {
        nndists.new[,s] = apply(nndist.sub[[s]], 1, sum)/length(ks) #average the ks distances
      }
      
      nndists = nndists.new
      
      datamarks = marks(datappp)
      uniquemarks = unique(datamarks)
      colnames(nndists) = unique(datamarks[datamarks!="Unknown"])
    }else{
      nndists  = nndists[,-which(colnames(nndists) == "Unknown")]
    }
    
	# Calculate the weights according the distances
    weight_num = apply(nndists, 1, min)/nndists
    init.weight = weight_num/apply(weight_num, 1, sum)
    init.weight[datamarks != "Unknown"] = rep(0, nclust)
    
	# Set the weights to 1 for the known points species and 0 when it is not the known points species
	# For the unknwon points the weights are the calculated weights
    for (i in 1:nclust)
    {
      rowfill = which(datamarks == colnames(nndists)[i])
      init.weight[rowfill, i] = 1
    }
    
    init.weight[is.nan(init.weight)] <- 1
    
    iterweights = init.weight
    itermarks = datamarks
    itermarks = colnames(nndists)[apply(init.weight, 1, which.max)]
  }
  
  #for kmeans method
  if (initweights == "kmeans"){
    # Get the coordinates of the clusters centers
	xy = coords(datappp)
    ncenter = nclust
    comp_mean = kmeans(xy, ncenter, nstart = nstart)
    
    Ccenter = comp_mean$centers 
    
	# calculate a distance of each point to each center
    All_dist_center = matrix(data=NA, nrow=nrow(xy), ncol=nclust)
    
    for (i in 1:(nclust))
    {
      Di= sqrt((xy$x-Ccenter[[i]])^2 + (xy$y-Ccenter[[i+nclust]])^2) 
      All_dist_center[,i] = Di
    }
    
    colnames(All_dist_center) <- paste("D", 1:nclust, sep = "")
    
    marksknown = uniquemarks[-which(uniquemarks == "Unknown")]
    colnames(All_dist_center) = marksknown
    
	# Calculate the weights according to the distances
    weight_num = apply(All_dist_center, 1, min)/All_dist_center
    init.weight = weight_num/apply(weight_num, 1, sum)
    
	# Set the weights to 1 for the known points species and 0 when it is not the known points species
	# For the unknwon points the weights are the calculated weights
    init.weight[datamarks != "Unknown"] = rep(0, nclust)
    for (i in 1:nclust)
    {
      rowfill = which(datamarks == colnames(All_dist_center)[i])
      init.weight[rowfill, i] = 1
    }
    
    init.weight[is.nan(init.weight)] <- 1
    
    iterweights = init.weight
    itermarks = datamarks
    itermarks = colnames(All_dist_center)[apply(init.weight, 1, which.max)]
    
  }
  
  # for random method
  if (initweights == "random"){
    BIC.mixt = c(rep(NA, nstart))
    init.weight = matrix(NA, datappp$n, nclust)
    
	# We randomly generate the weights of each observation for across the different species
    for (s in 1:nstart) {
      random.val = runif(nclust*datappp$n, min=0, max=1)
      weight_num = matrix(random.val, datappp$n, nclust)
      init.weight.rd = weight_num/apply(weight_num, 1, sum)  # make weights add up to 1
      
      colmarks = uniquemarks
      colmarks  = colmarks[-which(colmarks == "Unknown")]
      colnames(init.weight.rd) = colmarks
      
      # Set the weights to 1 for the known points species and 0 when it is not the known points species
	  # For the unknwon points the weights are the calculated weights
      init.weight.rd[datamarks != "Unknown"] = rep(0, nclust)
      for (i in 1:nclust)
      {
        rowfill = which(datamarks == colmarks[i])
        init.weight.rd[rowfill, i] = 1
      }
      init.weight.rd[is.nan(init.weight.rd)] <- 1
      
      iterweights = init.weight.rd
      itermarks = datamarks
      itermarks = colnames(init.weight.rd)[apply(init.weight.rd, 1, which.max)]
      
      # for random we fit a point process to choose through BIC crtierion the best random starting values
      p = table(itermarks)/sum(table(itermarks))
      
      iterppp = datappp
      marks(iterppp) = as.factor(itermarks)
      Qmask = makeMask(quads.)
      
      Q = quadscheme(data = iterppp, dummy = quads., method = "grid",
                     ntile = c(dim(Qmask)[2], dim(Qmask)[1]), 
                     npix = c(dim(Qmask)[2], dim(Qmask)[1]))
      
      formchr = as.character(ppmform)[2]
      formsplit = strsplit(formchr, "\\+")
      markform = as.formula(paste("~", paste(paste(formsplit[[1]], "* marks"), collapse = " + ")))
      
      if(is.null(kAreaInt)){
        fit1 = ppm(Q, trend = markform, covariates = cov.list., 
                   gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
      }else{
        fit1 = ppm(Q, trend = markform, covariates = cov.list., AreaInter(kAreaInt),
                   gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
        
      }
      
      # BIC criterion (Maasaro 2017 - Poisson mixture model selection)
      G = length(cov.list.)
      BIC.mixt[s]= -2*(fit1$maxlogpl) + (2*G-1)*(log(datappp$n) + 1)
      if(BIC.mixt[s] == min(BIC.mixt, na.rm = T)){
        init.weight = init.weight.rd
      }else{
        init.weight = init.weight
      }
      
      if(verbose) 
        cat(paste("Iteration", s, "\tBIC =", BIC.mixt[s], "\n"))
      
      s = s + 1
      
    }
    
    iterweights = init.weight
    itermarks = datamarks
    itermarks = colnames(init.weight)[apply(init.weight, 1, which.max)]
    
  }
  
  # random coin flip allocation system
  if (initweights == "CoinF"){
    itermarks = datamarks
    Lmarks  = length(datamarks[which(datamarks == "Unknown")])
    itermarks[which(datamarks == "Unknown")] = sample(c("Sp1", "Sp2", "Sp3"), Lmarks, replace=TRUE)
    
    colmarks = uniquemarks
    colmarks  = colmarks[-which(colmarks == "Unknown")]
    init.weight = matrix(0, nrow=length(itermarks) , ncol=nclust)
    for (i in 1:nclust)
    {
      rowfill = which(itermarks == colmarks[i])
      init.weight[rowfill, i] = 1
    }
    iterweights = init.weight
    colnames(iterweights)=c("Sp1", "Sp2", "Sp3")
  }
  
  #2# Fit point process models
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # continue general script for all the methods
  # initial proportion without unknown points
  p = table(datamarks[datamarks != "Unknown"])/sum(table(datamarks[datamarks != "Unknown"]))
  
  iterppp = datappp
  
  X=quads.$x
  Y=quads.$y
  quad.xy = data.frame(X, Y)
  
  ppp_list = n_known = keep_unk = keep_wt = list()
  Q = sp_wts.list = sp_list = list()
  n_unknown = splitppps[[length(splitppps)]]$n
  
  # We create some list of point pattern that have known points and all the unknown
  # Then we weight the quadrature point according to the known species (1 or 0) and the initial weights calculated
  for (i in 1:nclust) {
    n_known[[i]] = splitppps[[i]]$n
    keep_unk[[i]] = splitppps$Unknown[which(iterweights[unknown, i]!=0)]
    keep_wt[[i]] = iterweights[unknown, i][which(iterweights[unknown, i]!=0)]
    ppp_list[[i]] = superimpose(unmark(splitppps[[i]]), unmark(keep_unk[[i]]))
    
    Q[[i]]   = quadscheme(data = ppp_list[[i]], dummy = quads., method = "grid",
                          ntile = c(dim(Qmask)[2], dim(Qmask)[1]), npix = c(dim(Qmask)[2], dim(Qmask)[1]))
 
    sp_list[[i]] = data.frame(X = ppp_list[[i]]$x, Y=ppp_list[[i]]$y)

    sp_wts.list[[i]] = scoreweights(sp_list[[i]], quad.xy, 
                                    scores = c(rep(1, n_known[[i]]), keep_wt[[i]])) # generate quad weights from initial weights
  
    Q[[i]]$w = sp_wts.list[[i]]
    
  }
  
  # Fit Poisson PPMs
  fit1 = list()
  
  for (i in 1:nclust) {
    fit1[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list., 
                        gcontrol = list(epsilon = 1e-6, maxit = 100))
  }

  if(is.null(cov.bias)){
    cov.list. = cov.list.
  }else{
    pred.list = cov.list.
    set.Val = cov.bias #Variables to set to a certain value
    for (v in set.Val){
      pred.list[[v]]$v = kVal*pred.list[[v]]$v
    }
  }
  
  #continue script with bias
  formchr = as.character(ppmform)[2]
  #formsplit = strsplit(formchr, "\\+")
  #markform = as.formula(paste("~", paste(paste(formsplit[[1]], "* marks"), collapse = " + ")))
  
  #3# Compute predicted intensities
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fitbef.pred = pfit.b = list()
  for (i in 1:nclust) {
    if(is.null(cov.bias))
    {
      fitbef.pred[[i]] = predict(fit1[[i]], covariates = cov.list.)
    }
    else
    {
      fitbef.pred[[i]] = predict(fit1[[i]], covariates = pred.list)
    }
    
    if (plots == TRUE)
    {
      plot(fitbef.pred[[i]], main="predict - log(fitbef.pred)")
    }
    
    pfit.b[[i]] = fitted(fit1[[i]])
  }

  loglik.old <- 0.
  loglik.new <- 1.
  
  #
  # Iterator starts here, 
  #
  
  is_known = which(datamarks != "Unknown")
  niter <- 0
  while(abs(loglik.new - loglik.old)/abs(loglik.old) > tol) {
    if(niter >= maxit) {
      warning(paste("E-M algorithm failed to converge in",
                    maxit, ngettext(maxit, "iteration", "iterations")),
              call.=FALSE)
      break
    }
    
    
    niter <- niter + 1
    
    #4# Get the predicted intensities at the location S
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # E - step
    predint = matrix(NA, sum(unknown), nclust)
    for (i in 1:nclust)
    {
      ppp_i = ppp_Unknown
      #marks(ppp_i) = as.factor(colnames(iterweights)[i])
      if(is.null(cov.bias)){
        predint[,i] = predict(fit1[[i]], locations = ppp_i)
      }else{
        predint[,i] = predict(fit1[[i]], covariates = pred.list, locations = ppp_i)
      }
    }
    
    p_mat = t(matrix(p, ncol(predint), nrow(predint)))
    predint_p = predint*p_mat
    
    
    #5# Caluclate New weights
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    iterweights[unknown,] = predint_p/apply(predint_p, 1, sum)
    
    # end E-step
    
    #6# We compute species proportions
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # M - step
    p = apply(iterweights, 2, sum)/sum(apply(iterweights, 2, sum))
    
    itermarks = datamarks
    
    ppp_list_it = keep_unk_it = keep_wt_it = list()
    sp_wts_it.list = sp_list_it = list()

    # We create some list of point pattern that have known points and all the unknown
    # Then we weight the quadrature point according to the known species (1 or 0) and the initial weights calculated
    for (i in 1:nclust) {
      keep_unk_it[[i]] = splitppps$Unknown[which(iterweights[unknown, i]!=0)]
      keep_wt_it[[i]] = iterweights[unknown, i][which(iterweights[unknown, i]!=0)]
      ppp_list_it[[i]] = superimpose(unmark(splitppps[[i]]), unmark(keep_unk_it[[i]]))
      
      Q[[i]]   = quadscheme(data = ppp_list_it[[i]], dummy = quads., method = "grid",
                            ntile = c(dim(Qmask)[2], dim(Qmask)[1]), npix = c(dim(Qmask)[2], dim(Qmask)[1]))
      
      sp_list_it[[i]] = data.frame(X = ppp_list_it[[i]]$x, Y=ppp_list_it[[i]]$y)
      
      sp_wts_it.list[[i]] = scoreweights(sp_list_it[[i]], quad.xy, 
                                      scores = c(rep(1, n_known[[i]]), keep_wt_it[[i]])) # generate quad weights from initial weights
      
      Q[[i]]$w = sp_wts_it.list[[i]]
      
    }
  

    #8# Fit new point process models
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fit1.after =list()
    for (i in 1:n.sp) {
      if(is.null(kAreaInt)){
        fit1.after[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list., 
                         gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
      }else{
        fit1.after[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list., AreaInter(kAreaInt),
                         gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
      }
    }
    
    
    fit1 = fit1.after
    
    fitaft.pred = pfit.af = list()
    for (i in 1:n.sp) {
      if(is.null(cov.bias)){
        fitaft.pred[[i]]  = predict(fit1.after[[i]], covariates = cov.list.)
      }else{
        fitaft.pred[[i]]  = predict(fit1.after[[i]], covariates = pred.list)
      }
      
      if (plots == TRUE)
      {
        plot(envelope(datappp))
        plot(fitaft.pred[[i]], main="predict - log(fitaft.pred)")
      }
      
      pfit.af[[i]] = fitted(fit1.after[[i]])
    }

    #fitted.mix = fit1.after$internal$glmfit$fitted.values
    
    
    # end M-step
    
    #9# Stopping criterion
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # evaluate marginal loglikelihood
    
    loglik.old = loglik.new
    
    allp_mat = t(matrix(p, ncol(predint), nrow(iterweights)))
    loglik.new <- sum(log(apply(allp_mat * iterweights, 1, sum))) 
    
    loglik.new
    #loglik for mixture model : iterwights regroups weights for known and unknown observation
    
    # Weights for reclassified points 
    Newpts_w = iterweights[unknown,]
    
    sp_aug_ppp.list =  ppp_list
    
    # Predictions at locations
    #--
    if(is.null(cov.bias)){
      pred.loc = predict(fit1.after, locations = sp_int_im)
    }else{
      pred.loc = predict(fit1.after, covariates = pred.list, 
                         locations = sp_int_im)
    }
                           
    # Prepare weights for the plots
    Weight.df = as.data.frame(iterweights[unknown,])
    
    if(verbose) 
      cat(paste("Iteration", niter, "\tlogLik =", loglik.new,
                "\tp =", signif(p,4), "\n"))
  }
  
  if(verbose) {
    cat("\nEstimated parameters:\n")
    cat(paste("p [cluster] =", signif(p, 5), "\n"))
    cat(paste("\nloglik.new:\n", signif(loglik.new), "\n"))
    
    if (plots == TRUE)
    {
      par(xpd=NA)
      known.marks = unique(iterppp$marks)
      plot(x=seq_along(Weight.df[,1]), y=Weight.df[,1], col = "orange", pch=16, ylim=c(0,1),
           xlab="observations", ylab="weight")
      colvect=c("purple", "turquoise3", "darkred", "green", "brown")[1:nclust-1]
      for (i in 2:nclust) {
        points(x=seq_along(Weight.df[,i]), y=Weight.df[,i], col = colvect[i-1], pch=16, ylim=c(0,1))
        i =i + 1
        legend(200,1, c(known.marks), col = c("orange", colvect),
               pch = 16, xjust = 1, yjust = 0, merge = FALSE)
        
      }
      
    }
  }
  
  return(list(z = round(p, digits = 4),
              probs = p,
              niter = niter, maxit = maxit,
              converged = (niter >= maxit),
              New_weights = round(iterweights, digits = 4),
              pfit.b = pfit.b,
              pfit.af = pfit.af,
              #fitted.mix = fitted.mix,
              fit.final = fit1.after,
              fitaft.pred = fitaft.pred,
              pred.loc = pred.loc,
              Newpts_w = Newpts_w,
              iterppp = iterppp,
              sp_aug.list = sp_aug_ppp.list
              #hist=if(plothist) H else NULL
              # plot(x=seq_along(Weight.df$Sp1), y=Weight.df$Sp1, col = "orange", pch=16, ylim=c(0,1),
              #      xlab="observations", ylab="weight"),
              # points(x=seq_along(Weight.df$Sp2), y=Weight.df$Sp2, col = "purple", pch=18, ylim=c(0,1)),
              # points(x=seq_along(Weight.df$Sp3), y=Weight.df$Sp3, col = "Turquoise3", pch=17, ylim=c(0,1)),
              # legend(1,1, c("sp1", "sp2", "sp3"), col = c("orange", "purple", "Turquoise3"),
              #        pch = c(16, 18, 17), xjust = 1, yjust = 0, merge = FALSE)
  ))
}




### Function ppmLoopEngine
#---------------------------------------------------------------------

scoreweights = function(sp.xy, quad.xy, coord = c("X", "Y"), scores = NULL)
{
  if (is.null(scores)){
    score.all = rep(1, (dim(sp.xy)[1]) + dim(quad.xy)[1])
  }else{
    score.all = c(scores, rep(1, dim(quad.xy)[1]))
  }
  
  sp.col   = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) == coord[2]))
  quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy) == coord[2]))
  
  X.inc   = sort(unique(quad.xy[,quad.col[1]]))[2] - sort(unique(quad.xy[,quad.col[1]]))[1]
  Y.inc   = sort(unique(quad.xy[,quad.col[2]]))[2] - sort(unique(quad.xy[,quad.col[2]]))[1]
  quad.0X = min(quad.xy[,quad.col[1]]) - floor(min(quad.xy[,quad.col[1]])/X.inc)*X.inc
  quad.0Y = min(quad.xy[,quad.col[2]]) - floor(min(quad.xy[,quad.col[2]])/Y.inc)*Y.inc
  
  X = c(sp.xy[,quad.col[1]], quad.xy[,quad.col[1]])
  Y = c(sp.xy[,quad.col[2]], quad.xy[,quad.col[2]])
  
  round.X     = round((X - quad.0X)/X.inc)*X.inc
  round.Y     = round((Y - quad.0Y)/Y.inc)*Y.inc
  round.id    = paste(round.X, round.Y)
  round.tab   = aggregate(data.frame(score.all), list(ID = round.id), sum)
  scorewt     = X.inc*Y.inc*score.all/round.tab$score.all[match(round.id, round.tab$ID)]
  scorewt
}

#------------------------------------------------------------------------------------
#  								function ppmAddEngine
#------------------------------------------------------------------------------------

ppmLoopEngine = function(datappp, all_test, n.sp, addpt = c("LoopA","LoopT", "LoopE"), quads., win.,
                         ppmform, delta_max=NULL, delta_min=NULL, delta_step =NULL, num.add = NULL,
                         cov.list., cov.bias=NULL, kVal =NULL, kAreaInt=NULL, maxit = 50,
                         tol=0.000001, verbose = TRUE, plots = FALSE){
  
  # Define some objects
  datamarks = marks(datappp)
  uniquemarks = unique(datamarks)
  unknown = datamarks == "Unknown"
  nclust  = length(unique(datamarks)) - 1
  Qmask = makeMask(quads.)
  
  # Separate the multi type point pattern
  splitppps = split(datappp, as.factor(marks(datappp)))
  for (i in 1:(nclust + 1))
  {
    assign(paste("sp_sub", names(splitppps)[i], sep = ""),
           splitppps[[i]])
  }
  
  #1# Fit initial point processes
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Specie fitted separetely
  
  ppp_list = list()
  Q = list()
  
  for (i in 1:nclust) {
    ppp_list[[i]] = unmark(splitppps[[i]])
    Q[[i]]   = quadscheme(data = ppp_list[[i]], dummy = quads., method = "grid",
                          ntile = c(dim(Qmask)[2], dim(Qmask)[1]), npix = c(dim(Qmask)[2], dim(Qmask)[1]))
    i=i+1
  }
  
  # Fit Poisson PPMs
  ppm_list = list()
    
  for (i in 1:nclust) {
    ppm_list[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list, 
                        gcontrol = list(epsilon = 1e-6, maxit = 100))
    i=i+1
  }
  
  X=quads.$x
  Y=quads.$y
  quad.xy = data.frame(X, Y)
  
  if(is.null(cov.bias)){
    cov.list. = cov.list.
  }else{#--- Set observer bias variables to kVal 
    pred.list = cov.list.
    set.Val = cov.bias #Variables to set to a certain value
    for (v in set.Val){
      pred.list[[v]]$v = kVal*pred.list[[v]]$v
    }
  }
  
  datamarks = marks(datappp)
  uniquemarks = unique(datamarks)
  unknown = datamarks == "Unknown"
  names.mark = uniquemarks[uniquemarks != "Unknown"]
  
  niter <- 0
  
  loglik.old.sp = rep(NA, nclust)
  loglik.new.sp = rep(NA, nclust)
  
  for (i in 1:nclust) {
    loglik.old.sp[i] <- 1.
    loglik.new.sp[i] = ppm_list[[i]]$maxlogpl
  }
  
  Lcrit = 1.
  Lcrit.vec = rep(NA, maxit)
  breakloop = 0
  
  is_known = which(datamarks != "Unknown")
  
  all_wts = array(data = NA, dim = c((maxit + 1), all_test$n, n.sp))
  
  while(breakloop == 0)
  {
    niter = niter + 1
	
    #2 Compute predicted intensities
    
    pr_ppm_list = list()
    for (i in 1:nclust) {
      if(is.null(cov.bias))
      {
        pr_ppm_list[[i]] = predict(ppm_list[[i]], locations = all_test)
      }
      else
      {
        pr_ppm_list[[i]] = predict(ppm_list[[i]], covariates = pred.list, locations = all_test)
      }
    }
    
    
    #3 Compute membership probabilities
    
    all_preds = data.frame(matrix(unlist(pr_ppm_list),
                                  nrow=length(pr_ppm_list[[1]]), byrow=F))
    
    test_wts  = all_preds/apply(all_preds, 1, sum)
    all_wts[niter,,] = as.matrix(test_wts)
    
    max_pred = apply(all_preds, 1, which.max)
    
    max_pred.vec = rep(NA, nclust)
    for (i in 1:nclust) {
      max_pred.vec[i] = sum(max_pred == i)
      i=i+1
    }

    pred.check = as.vector(max_pred.vec)
    # Set up otpion for initial weights  
    addpt <- match.arg(addpt)
    
    #4 Augment points for point patterns
	
    if (addpt == "LoopA")
    {
      addtosp.list =list()
      for (i in 1:nclust) {
        addtosp.list[[i]] = (1:all_test$n)
        i=i+1
      }

    }
    if (addpt == "LoopT")
    {
      addtosp.list =list()
      for (i in 1:nclust) {
        addtosp.list[[i]] = which(test_wts[,i] > delta_max)
        i=i+1
      }
      
    }
    if (addpt == "LoopE")
    {
      if(num.add > all_test$n/n.sp)
      {
        print("Impossible to add so many points, the highest possible number will be used instead")
        num.add = floor(all_test$n/n.sp)
      }
      else
      {
        num.add = num.add
      }
      
      add_max = apply(test_wts, 2, sort, decreasing = TRUE)[num.add,]
      addtosp.list =list()
      for (i in 1:nclust) {
        addtosp.list[[i]] = if(anyNA(all_test$x[test_wts[,i] >= add_max[i]]) == TRUE) integer() else which(test_wts[,i] >= add_max[i])
        i=i+1
      }
      
    }
    
    # lists and vectors needed in the next steps
    sp_aug.list = sp_wts.list = sp_aug_ppp.list = Q_aug.list = ppm.L.list = list()
    ppm.L.pred = ppm.pred.list = list()
    Dloglik = counts.sp = rep(NA, nclust)
    
    for (i in 1:nclust) {
      sp_aug.list[[i]] = data.frame(X = c(ppp_list[[i]]$x, all_test$x[addtosp.list[[i]]]), Y = c(ppp_list[[i]]$y, all_test$y[addtosp.list[[i]]])) # add unknown points to known points of species 1
      quad.xy = data.frame(X, Y)
      
      #5 update quadrature weights
      
      #scores for species with known label (weight =1) and for the new obs (test_wts)
      sp_wts.list[[i]] = scoreweights(sp_aug.list[[i]], quad.xy, scores = c(rep(1, ppp_list[[i]]$n), test_wts[addtosp.list[[i]], i])) # generate quad weights for augmented species 1
      
      # Augmented point patterns
      win.   = win. # added because HPC wouldn't work
      
      sp_aug_ppp.list[[i]] = ppp(x =sp_aug.list[[i]]$X, y = sp_aug.list[[i]]$Y, window = win.)
      
      # Augmented quadrature scheme
      Q_aug.list[[i]] = quadscheme(data = sp_aug_ppp.list[[i]], dummy = quads., method = "grid",
                                   ntile = c(dim(Qmask)[2], dim(Qmask)[1]), npix = c(dim(Qmask)[2], dim(Qmask)[1]))
      
      # Replace quadrature weights with those calculated with the scoreweights function
      # This is necessary because spatstat's quadscheme function treats all points the same.
      # We want to treat the points with unknown labels as "fractional" points with weights coming from the single-species PPMs
      
      Q_aug.list[[i]]$w = sp_wts.list[[i]]
      
      #6# Fit new point processes using the augmented points patterns and the quadrature weights
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # Augmented PPMs
      
      if(is.null(kAreaInt))
      {
        ppm.L.list[[i]] = ppm(Q_aug.list[[i]], trend = ppmform, covariates = cov.list., 
                              gcontrol = list(epsilon = 1e-6, maxit = 100))
      }
      else
      {
        ppm.L.list[[i]] = ppm(Q_aug.list[[i]], trend = ppmform, covariates = cov.list.,
                              AreaInter(kAreaInt),
                              gcontrol = list(epsilon = 1e-6, maxit = 100))
      }
      
      if (plots == TRUE)
      {
        if(is.null(cov.bias))
        {
          ppm.L.pred[[i]]  = predict(ppm.L.list[[i]])
        }
        else
        {
          ppm.L.pred[[i]]  = predict(ppm.L.list[[i]], covariates = pred.list)

        }
        
        plot(ppm.L.pred[[i]], main="predict - ppm.pred") 
      }
      
      # to get the weights
      if(is.null(cov.bias))
      {
        ppm.pred.list[[i]]  = predict(ppm.L.list[[i]], location=datappp)
      }
      else
      {
        ppm.pred.list[[i]]  = predict(ppm.L.list[[i]], covariates = pred.list, location=datappp)
      }
      
      # counts per species
      counts.sp[i] = ppp_list[[i]]$n + sum(test_wts[,i])
      
      #7 Stopping criterion
      
      loglik.old.sp[i] = loglik.new.sp[i]
      loglik.new.sp[i] = ppm.L.list[[i]]$maxlogpl

      Dloglik[i] = abs(loglik.new.sp[i] - loglik.old.sp[i])
      i=i+1
    }
    
    DiffL = sum(Dloglik)
    sumL.new = abs(sum(loglik.new.sp))
    
    Lcrit = DiffL/sumL.new 
    Lcrit.vec[niter] = Lcrit

    itercounts = counts.sp
    p = itercounts/sum(itercounts)
    
    if(verbose)
    {
      cat(paste("Iteration", niter, "\tLcrit =", Lcrit,
                "\tp =", signif(p,4), "\n"))
    }
    
    # redefine ppms for next iteration
    for (i in 1:nclust) {
      ppm_list[[i]] = ppm.L.list[[i]]
    }
    
    # break loop
    
    if (Lcrit < tol)
    {
      breakloop = 1
    }
    
    if (niter == maxit)
    {
      breakloop = 1
    }
    
    if (addpt == "LoopT")
    {
      delta_max = delta_max - delta_step
      if (delta_max < delta_min)
	  {
		breakloop = 1
      }
    }
    
    if (addpt == "LoopE")
    {
      num.add = num.add + 1
      if (num.add > all_test$n/n.sp)
      {
        breakloop = 1
      }
    }
    
  }
  
  # Compute final predicted intensities
  pr_ppm_list.unk = list()
  for (i in 1:nclust) {
    if(is.null(cov.bias))
    {
      pr_ppm_list[[i]] = predict(ppm_list[[i]], locations = datappp)
    }
    else
    {
      pr_ppm_list[[i]] = predict(ppm_list[[i]], covariates = pred.list, locations = datappp)
    }
    
    pr_ppm_list.unk[[i]] = pr_ppm_list[[i]][-is_known]
  }
  
  #3 Compute final membership probabilities
  
  all_preds = data.frame(matrix(unlist(pr_ppm_list.unk),
                                nrow=length(pr_ppm_list.unk[[1]]), byrow=F))
  
  test_wts  = all_preds/apply(all_preds, 1, sum)
  all_wts[niter + 1,,] = as.matrix(test_wts)
  
  
  # prediction at locations
  #--
  pred.loc = list()
  for (l in 1:n.sp) {
    if(is.null(cov.bias)){
      pred.loc[[l]] = predict(ppm_list[[l]], locations = quads.win)
    }else{
      pred.loc[[l]] = predict(ppm_list[[l]], covariates = pred.list, 
                                        locations = quads.win)
    }
  }
  
  if (plots == TRUE)
  {
    par(xpd=NA)
    known.marks = unique(iterppp$marks)
    plot(x=seq_along(Weight.df[,1]), y=Weight.df[,1], col = "orange", pch=16, ylim=c(0,1),
         xlab="observations", ylab="weight")
    colvect=c("purple", "turquoise3", "darkred", "green", "brown")[1:nclust-1]
    for (i in 2:nclust) {
      points(x=seq_along(Weight.df[,i]), y=Weight.df[,i], col = colvect[i-1], pch=16, ylim=c(0,1))
      i =i + 1
      legend(110,1, c(known.marks), col = c("orange", colvect),
             pch = 16, xjust = 1, yjust = 0, merge = FALSE)
      
    }
    
  }
  
  return(list(z = round(p, digits = 4),
              New_weights = test_wts,
              Newpts_w = test_wts,
              all_wts = all_wts,
              ppm_list = ppm_list,
              niter = niter, 
              ppm.pred.list = pr_ppm_list,
              ppm.pred.list.unk = pr_ppm_list.unk,
              sp_aug.list = sp_aug.list,
              sp_aug_ppp.list = sp_aug_ppp.list,
              delta_max = delta_max,
              delta_min = delta_min,
              delta_step = delta_step,
              pred.loc = pred.loc,
              sp_aug.list = sp_aug.list
              #hist=if(plothist) H else NULL
              # plot(x=seq_along(Weight.df$Sp1), y=Weight.df$Sp1, col = "orange", pch=16, ylim=c(0,1),
              #      xlab="observations", ylab="weight"),
              # points(x=seq_along(Weight.df$Sp2), y=Weight.df$Sp2, col = "purple", pch=18, ylim=c(0,1)),
              # points(x=seq_along(Weight.df$Sp3), y=Weight.df$Sp3, col = "Turquoise3", pch=17, ylim=c(0,1)),
              # legend(1,1, c("sp1", "sp2", "sp3"), col = c("orange", "purple", "Turquoise3"),
              #        pch = c(16, 18, 17), xjust = 1, yjust = 0, merge = FALSE)
  ))
}



##------------------------------------------------------------------------------
#              functions and measures of performance
#------------------------------------------------------------------------------
###----------------- IMSE

IMSE = function(mu1, mu2, fun = "log", mu.min = 1.e-5, rescale = TRUE)
{
  mu1.use = mu1
  mu1.use[mu1.use < mu.min] = mu.min
  mu2.use = mu2
  mu2.use[mu2.use < mu.min] = mu.min
  
  if (rescale == TRUE)
  {
    mu2.use = mu2.use*mean(mu1.use)/mean(mu2.use)
  }
  
  if (fun == "log")
  {
    mu1.use = log(mu1.use)
    mu2.use = log(mu2.use)
  }
  if (fun == "sqrt")
  {
    mu1.use = sqrt(mu1.use)
    mu2.use = sqrt(mu2.use)
  }
  imse = sum((mu1.use - mu2.use)^2)
  imse
}

###----------------- corint for sumcor calculation

corint = function(mu1, mu2, fun = "log", method=c("pearson", "kendall", "spearman"), mu.min = 1.e-5)
{
  mu1.use = mu1
  mu1.use[mu1.use < mu.min] = mu.min
  mu2.use = mu2
  mu2.use[mu2.use < mu.min] = mu.min
  if (fun == "log")
  {
    mu1.use = log(mu1.use)
    mu2.use = log(mu2.use)
  }
  
  # Set up otpion for initial weights  
  addpt <- match.arg(method)
  if (method == "pearson"){
    corint.pea = cor(mu1.use, mu2.use, method = "pearson")
    return(corint.pea)
  }
  
  if(method == "kendall"){
    corint.kend = cor(mu1.use, mu2.use, method = "kendall")
    return(corint.kend)
  }
  
  if(method == "spearman"){
    corint.spea = cor(mu1.use, mu2.use, method = "spearman")
    return(corint.spea)
  }
  
}

###----------------- RSS

RSS = function(weightmatrix, test_labels)
{
  mark.cols = match(test_labels, colnames(weightmatrix))
  correctweights = weightmatrix[cbind(seq_along(mark.cols), mark.cols)]
  RSSscore = sum((correctweights - 1)^2)
  RSSscore
}

###----------------- Accuracy

Accuracy = function(all_true, New_weights, test_labels, n.sp){
  W.max = apply(New_weights[(1:nrow(New_weights)),], 1, max)
  C.id = apply(New_weights[(1:nrow(New_weights)),],
               1,function(x) which(x==max(x)))
  
  indiv_testlab = as.data.frame(cbind(test_labels, C.id, W.max))
  
  levels(indiv_testlab$test_labels) <- c(unique(test_labels))
  
  # New method for accuracy
  levels(indiv_testlab$test_labels)
  levels(indiv_testlab$C.id)
  for (i in 1:n.sp) {
    levels(indiv_testlab$test_labels)[levels(indiv_testlab$test_labels)==paste("Sp", i, sep = "")] <- paste(i, sep = "")
    
  }
  
  acc_tab <- table(factor(indiv_testlab$test_labels, levels=1:n.sp), 
                   factor(indiv_testlab$C.id, levels=1:n.sp))
  
  #dim(acc_tab)
  
  CM.acc = confusionMatrix(acc_tab)
  m.acc= CM.acc$table
  
  # some usuful calc
  n = sum(m.acc) # number of instances
  nc = nrow(m.acc) # number of classes
  diag = diag(m.acc) # number of correctly classified instances per class 
  rowsums = apply(m.acc, 1, sum) # number of instances per class
  colsums = apply(m.acc, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  
  # accuracy measure
  acc = sum(diag) / n 
  acc 
  
}


Perffunc = function(fit, sp_int.list, datappp, fun = "log", rescale = TRUE, 
                    method=c("pearson", "kendall", "spearman"), LoopM=TRUE,
                    mu.min = 1.e-5, all_true, test_labels, n.sp, pf = c(NULL)){
  
  pred = fit$pred.loc
  New_weights = fit$Newpts_w
  
  if(LoopM==TRUE){
    colnames(New_weights)= c(unique(test_labels))
  }
  
  ###----------------- IMSE
  
  IMSE = function(mu1, mu2, fun = "log", mu.min = 1.e-5, rescale = TRUE)
  {
    mu1.use = mu1
    mu1.use[mu1.use < mu.min] = mu.min
    mu2.use = mu2
    mu2.use[mu2.use < mu.min] = mu.min
    
    if (rescale == TRUE)
    {
      mu2.use = mu2.use*mean(mu1.use)/mean(mu2.use)
    }
    
    if (fun == "log")
    {
      mu1.use = log(mu1.use)
      mu2.use = log(mu2.use)
    }
    if (fun == "sqrt")
    {
      mu1.use = sqrt(mu1.use)
      mu2.use = sqrt(mu2.use)
    }
    imse = sum((mu1.use - mu2.use)^2)
    imse
  }
  
  ###----------------- corint for sumcor calculation
  
  corint = function(mu1, mu2, fun = "log", method=c("pearson", "kendall", "spearman"), mu.min = 1.e-5)
  {
    mu1.use = mu1
    mu1.use[mu1.use < mu.min] = mu.min
    mu2.use = mu2
    mu2.use[mu2.use < mu.min] = mu.min
    if (fun == "log")
    {
      mu1.use = log(mu1.use)
      mu2.use = log(mu2.use)
    }
    
    # Set up otpion for initial weights  
    addpt <- match.arg(method)
    if (method == "pearson"){
      corint.pea = cor(mu1.use, mu2.use, method = "pearson")
      return(corint.pea)
    }
    
    if(method == "kendall"){
      corint.kend = cor(mu1.use, mu2.use, method = "kendall")
      return(corint.kend)
    }
    
    if(method == "spearman"){
      corint.spea = cor(mu1.use, mu2.use, method = "spearman")
      return(corint.spea)
    }
    
  }
  
  ###----------------- RSS
  
  RSS = function(weightmatrix, test_labels)
  {
    mark.cols = match(test_labels, colnames(weightmatrix))
    correctweights = weightmatrix[cbind(seq_along(mark.cols), mark.cols)]
    RSSscore = sum((correctweights - 1)^2)
    RSSscore
  }
  
  ###----------------- Accuracy
  
  Accuracy = function(all_true, New_weights, test_labels, n.sp){
    W.max = apply(New_weights[(1:nrow(New_weights)),], 1, max)
    C.id = apply(New_weights[(1:nrow(New_weights)),],
                 1,function(x) which(x==max(x)))
    
    indiv_testlab = as.data.frame(cbind(test_labels, C.id, W.max))
    
    levels(indiv_testlab$test_labels) <- c(unique(test_labels))
    
    # New method for accuracy
    levels(indiv_testlab$test_labels)
    levels(indiv_testlab$C.id)
    for (i in 1:n.sp) {
      levels(indiv_testlab$test_labels)[levels(indiv_testlab$test_labels)==paste("Sp", i, sep = "")] <- paste(i, sep = "")
      
    }
    
    acc_tab <- table(factor(indiv_testlab$test_labels, levels=1:n.sp), 
                     factor(indiv_testlab$C.id, levels=1:n.sp))
    
    #dim(acc_tab)
    
    CM.acc = confusionMatrix(acc_tab)
    m.acc= CM.acc$table
    
    # some usuful calc
    n = sum(m.acc) # number of instances
    nc = nrow(m.acc) # number of classes
    diag = diag(m.acc) # number of correctly classified instances per class 
    rowsums = apply(m.acc, 1, sum) # number of instances per class
    colsums = apply(m.acc, 2, sum) # number of predictions per class
    p = rowsums / n # distribution of instances over the actual classes
    q = colsums / n # distribution of instances over the predicted classes
    
    # accuracy measure
    acc = sum(diag) / n 
    acc 
    
  }
  
  
  if(any(pf)=="IMSE"){
    sp.predlist = list()
    for (l in 1:n.sp) {

      if (LoopM == TRUE){
        sp.predlist[[l]] = as.vector(t(pred[[l]]))
        IMSEscore = sum(IMSE(sp_int.list[[l]], (sp.predlist[[l]]/(length(fit$sp_aug.list[[l]]$X)))*datappp$n/n.sp))
        
      }else{
        sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
        IMSEscore = sum(IMSE(sp_int.list[[l]], (sp.predlist[[l]]/(fit$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        
      }
    }
  }
  if(any(pf)=="corint"){
    sp.predlist = list()
    for (l in 1:n.sp) {
      sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
      
      sumcor1 = sum(corint(sp_int.list[[l]], sp.predlist[[l]], method="pearson"))
      sumcor2 = sum(corint(sp_int.list[[l]], sp.predlist[[l]], method="kendall"))
      
    }
  }
  
  if(any(pf)=="mRSS"){
    RSS = RSS(New_weights, test_labels)
    meanRSS = RSS(New_weights, test_labels)/length(test_labels)
  }
  
  if(any(pf)=="Acc"){
    accmat = Accuracy(all_true, New_weights, test_labels, n.sp)
  }
  
  if(is.null(pf)==TRUE){
    accmat = Accuracy(all_true, New_weights, test_labels, n.sp)
    RSSscore = RSS(New_weights, test_labels)
    meanRSS = RSS(New_weights, test_labels)/length(test_labels)
    
    sp.predlist = list()
    for (l in 1:n.sp) {
      if (LoopM == TRUE){
        sp.predlist[[l]] = as.vector(t(pred[[l]]))
        IMSEscore = sum(IMSE(sp_int.list[[l]], (sp.predlist[[l]]/(length(fit$sp_aug.list[[l]]$X)))*datappp$n/n.sp))
        
      }else{
        sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
        IMSEscore = sum(IMSE(sp_int.list[[l]], (sp.predlist[[l]]/(fit$sp_aug.list[[l]]$n))*datappp$n/n.sp))
        
      }
      
      sumcor1 = sum(corint(sp_int.list[[l]], sp.predlist[[l]], method="pearson"))
      sumcor2 = sum(corint(sp_int.list[[l]], sp.predlist[[l]], method="kendall"))
      
    }
  }
  
  return(list(accmat=accmat, RSS=RSSscore, meanRSS=meanRSS, IMSE=IMSEscore, sumcor1=sumcor1, sumcor2=sumcor2))
}


Predlist = function(pred, sp){
  sp.predlist = list()
  for (l in 1:n.sp) {
    sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
  }
  return(sp.predlist)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to combine the different methods to compare
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
