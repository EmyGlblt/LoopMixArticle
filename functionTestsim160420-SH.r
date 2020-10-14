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

ppmMixEngine = function(Known.ppp, Unknown.ppp, quadsenv, ppmform,  
                        initweights = c("knn", "kps", "kmeans", "random", "CoinF"), cov.list,
                        k=1, ks=1, nstart=nstart, cov.bias=NULL, SetValBias = NULL, rAreaInter=NULL,
                        verbose = TRUE, tol = 0.00001, maxit = 50, plots = FALSE, classif = c("hard", "soft"))
{
  
  n.sp = length(unique(marks(Known.ppp)))
  
  datappp = superimpose.ppp(Known.ppp, Unknown.ppp)
  # define some objects
  datamarks = marks(datappp)
  uniquemarks = unique(datamarks)
  unknown = datamarks == "Unknown"
  names.mark = uniquemarks[uniquemarks != "Unknown"]
  nclust  = length(unique(datamarks)) - 1
  
  # prepare quadrature ppp
  winMix   = owin(xrange = c(min(quadsenv$X)-0.5, max(quadsenv$X)+0.5), 
               yrange = c(min(quadsenv$Y)-0.5, max(quadsenv$Y)+0.5))
  quads. = ppp(quadsenv$X, quadsenv$Y, window = winMix)
  Qmask = makeMask(quads.)
  
  # prepare list of images of covariates
  cov.list = list()
  names.env = c()
  for (v in 1:(length(quadsenv)-2))
  {
    v.v = as.im(data.frame(x = quadsenv$X, y = quadsenv$Y, z = quadsenv[,v+2]))
    cov.list[[v]] = v.v
    names.env[v] = names(quadsenv)[v+2]
  }
  names(cov.list) = names.env
  
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
    if(k > Known.ppp$n)
    {
      print("Impossible to consider so many points, the highest possible number will be used instead")
      k = Known.ppp$n
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
    if(max(ks) > min(table(Known.ppp$marks)))
    {
      print("Impossible to consider so many points, the highest possible number will be used instead")
      ks = 1:min(table(Known.ppp$marks))
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
    init.weight = weight_num/apply(weight_num, 1, sum)     #IR: This will give all points (both with known and unknown species labels) weights between 0 and 1
    
	# Set the weights to 1 for the known points species and 0 when it is not the known points species
	# For the unknown points the weights are the calculated weights
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
    #BIC.mixt = c(rep(NA, nstart))
    init.weight = matrix(NA, datappp$n, nclust)
    
	# We randomly generate the weights of each observation for across the different species
    #for (s in 1:nstart) {
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
      
      if(is.null(rAreaInter)){
        fit1 = ppm(Q, trend = markform, covariates = cov.list, 
                   gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
      }else{
        fit1 = ppm(Q, trend = markform, covariates = cov.list, AreaInter(rAreaInter),
                   gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
        
      }
      
      # # BIC criterion (Maasaro 2017 - Poisson mixture model selection)
      # G = length(cov.list)
      # BIC.mixt[s]= -2*(fit1$maxlogpl) + (2*G-1)*(log(datappp$n) + 1)
      # if(BIC.mixt[s] == min(BIC.mixt, na.rm = T)){
        # init.weight = init.weight.rd
      # }else{
        # init.weight = init.weight
      # }
      
      # if(verbose) 
        # cat(paste("Iteration", s, "\tBIC =", BIC.mixt[s], "\n"))
      
      # s = s + 1
      
    #}
    
    iterweights = init.weight.rd
    itermarks = datamarks
    itermarks = colnames(init.weight.rd)[apply(init.weight.rd, 1, which.max)]
    
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
  
  
  
  ## hard classification
  # Set up otpion for initial weights  
  classif <- match.arg(classif)
  
  if(classif == "hard"){
    #2# Fit point process models
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # continue general script for all the methods
  p = table(itermarks)/sum(table(itermarks))
  
  iterppp = datappp
  marks(iterppp) = as.factor(itermarks)
  Qmask = makeMask(quads.)
  
  Q = quadscheme(data = iterppp, dummy = quads., method = "grid",
                 ntile = c(dim(Qmask)[2], dim(Qmask)[1]), 
                 npix = c(dim(Qmask)[2], dim(Qmask)[1]))
  
  if(is.null(cov.bias)){
    cov.list = cov.list
  }else{
    pred.list = cov.list
    set.Val = cov.bias #Variables to set to a certain value
    for (v in set.Val){
      pred.list[[v]]$v = SetValBias*pred.list[[v]]$v
    }
  }
  
  #continue script with bias
  formchr = as.character(ppmform)[2]
  formsplit = strsplit(formchr, "\\+")
  markform = as.formula(paste("~", paste(paste(formsplit[[1]], "* marks"), collapse = " + ")))
  
  if(is.null(rAreaInter)){
    fit1 = ppm(Q, trend = markform, covariates = cov.list, 
               gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
  }else{
    fit1 = ppm(Q, trend = markform, covariates = cov.list, AreaInter(rAreaInter),
               gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
    
  }
  
  #3# Compute predicted intensities
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  if(is.null(cov.bias)){
    fitbef.pred  = predict(fit1, covariates = cov.list)
  }else{
    fitbef.pred  = predict(fit1, covariates = pred.list)
  }
  
  if (plots == TRUE)
  {
    plot(fitbef.pred, main="predict - log(fitbef.pred)")
  }
  
  pfit.b = fitted(fit1)
  
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
      marks(ppp_i) = as.factor(colnames(iterweights)[i])
      if(is.null(cov.bias)){
        predint[,i] = predict(fit1, locations = ppp_i)
      }else{
        predint[,i] = predict(fit1, covariates = pred.list, locations = ppp_i)
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
    itermarks = colnames(iterweights)[apply(iterweights, 1, which.max)] # assign marks based on new weights
    
    iterppp = datappp
    marks(iterppp) = as.factor(itermarks)
    
    
    #7# Update quadrature weights
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Q = quadscheme(data = iterppp, dummy = quads., method = "grid",
                   ntile = c(dim(Qmask)[2], dim(Qmask)[1]),
                   npix = c(dim(Qmask)[2], dim(Qmask)[1]))
    
    formchr = as.character(ppmform)[2]
    formsplit = strsplit(formchr, "\\+")
    markform = as.formula(paste("~", paste(paste(formsplit[[1]], "* marks"), collapse = " + ")))
    
    #Q$w = sp_wts
    
    #8# Fit new point process models
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if(is.null(rAreaInter)){
      fit1.after = ppm(Q, trend = markform, covariates = cov.list, 
                       gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
    }else{
      fit1.after = ppm(Q, trend = markform, covariates = cov.list, AreaInter(rAreaInter),
                       gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
    }
    
    fit1 = fit1.after
    
    if(is.null(cov.bias)){
      fitaft.pred  = predict(fit1.after, covariates = cov.list)
    }else{
      fitaft.pred  = predict(fit1.after, covariates = pred.list)
    }
    
    
    if (plots == TRUE)
    {
      #plot(envelope(iterppp))
      plot(fitaft.pred, main="predict - log(fitaft.pred)")
    }
    
    pfit.af = fitted(fit1.after)
    
    fitted.mix = fit1.after$internal$glmfit$fitted.values
    
    #m.cor <- markcorr(iterppp)
    
    #if (plots == TRUE)
    #{
    #  plot(m.cor)
    #}
    
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
    
    sp_aug_ppp.list = iterppp

    
    # Predictions at locations
    #--
    if(is.null(cov.bias)){
      pred.loc = predict(fit1.after, locations = cov.list[[1]])
    }else{
      pred.loc = predict(fit1.after, covariates = pred.list,
                         locations = cov.list[[1]])
    }
    
    
    # Prepare weights for the plots
    Weight.df = as.data.frame(iterweights[unknown,])
    
    
    if(verbose) 
      cat(paste("Iteration", niter, "\tlogLik =", loglik.new,
                "\tp =", signif(p,4), "\n"))
  }
  }
  
  if(classif == "soft"){
  #2# Fit point process models
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # continue general script for all the methods
  # initial proportion without unknown points
  p = apply(iterweights, 2, sum)/sum(apply(iterweights, 2, sum))  
  
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
    fit1[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list, 
                        gcontrol = list(epsilon = 1e-6, maxit = 100))
  }

  if(is.null(cov.bias)){
    cov.list = cov.list
  }else{
    pred.list = cov.list
    set.Val = cov.bias #Variables to set to a certain value
    for (v in set.Val){
      pred.list[[v]]$v = SetValBias*pred.list[[v]]$v
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
      fitbef.pred[[i]] = predict(fit1[[i]], covariates = cov.list)
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
      if(is.null(rAreaInter)){
        fit1.after[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list, 
                         gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
      }else{
        fit1.after[[i]] = ppm(Q[[i]], trend = ppmform, covariates = cov.list, AreaInter(rAreaInter),
                         gcontrol = list(epsilon = 1e-6, maxit = 100)) # including known and unknown points
      }
    }
    
    
    fit1 = fit1.after
    
    fitaft.pred = pfit.af = list()
    for (i in 1:n.sp) {
      if(is.null(cov.bias)){
        fitaft.pred[[i]]  = predict(fit1.after[[i]], covariates = cov.list)
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
    pred.loc=list()
    for (i in 1:n.sp) {
      if(is.null(cov.bias)){
        pred.loc[[i]]  = predict(fit1.after[[i]], locations = cov.list[[1]])
      }else{
        pred.loc[[i]]  = predict(fit1.after[[i]], covariates = pred.list,
                                 locations = cov.list[[1]])
      }
    }
    
                           
    # Prepare weights for the plots
    Weight.df = as.data.frame(iterweights[unknown,])
    
    if(verbose) 
      cat(paste("Iteration", niter, "\tlogLik =", loglik.new,
                "\tp =", signif(p,4), "\n"))
  }
  
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
              classif = classif,
              New_weights = round(iterweights, digits = 4),
              Weight.df = Weight.df,
              pfit.b = pfit.b,
              pfit.af = pfit.af,
              #fitted.mix = fitted.mix,
              fit.final = fit1.after,
              fitaft.pred = fitaft.pred,
              pred.loc = pred.loc,
              Newpts_w = Newpts_w,
              iterppp = iterppp,
              sp_aug.list = sp_aug_ppp.list,
              Known.ppp = Known.ppp,
              Unknown.ppp =Unknown.ppp
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

ppmLoopEngine = function(Known.ppp, Unknown.ppp, addpt = c("LoopA","LoopT", "LoopE"), quadsenv,
                         ppmform, delta_max=0.9, delta_min=0.5, delta_step =0.1, num.add = 1,
                         cov.bias=NULL, SetValBias =NULL, rAreaInter=NULL, maxit = 50,
                         tol=0.000001, verbose = TRUE, plots = FALSE){
  
  n.sp = length(unique(marks(Known.ppp)))
  
  datappp = superimpose.ppp(Known.ppp, Unknown.ppp)
  # Define some objects
  datamarks = marks(datappp)
  uniquemarks = unique(datamarks)
  unknown = datamarks == "Unknown"
  nclust  = length(unique(datamarks)) - 1

  # Separate the multi type point pattern
  splitppps = split(datappp, as.factor(marks(datappp)))
  for (i in 1:(nclust + 1))
  {
    assign(paste("sp_sub", names(splitppps)[i], sep = ""),
           splitppps[[i]])
  }
  
  # prepare quadrature ppp
  winMix   = owin(xrange = c(min(quadsenv$X)-0.5, max(quadsenv$X)+0.5), 
                  yrange = c(min(quadsenv$Y)-0.5, max(quadsenv$Y)+0.5))
  quads. = ppp(quadsenv$X, quadsenv$Y, window = winMix)
  Qmask = makeMask(quads.)
  
  # prepare list of images of covariates
  cov.list = list()
  names.env = c()
  for (v in 1:(length(quadsenv)-2))
  {
    v.v = as.im(data.frame(x = quadsenv$X, y = quadsenv$Y, z = quadsenv[,v+2]))
    cov.list[[v]] = v.v
    names.env[v] = names(quadsenv)[v+2]
  }
  names(cov.list) = names.env
  
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
    cov.list = cov.list
  }else{#--- Set observer bias variables to SetValBias 
    pred.list = cov.list
    set.Val = cov.bias #Variables to set to a certain value
    for (v in set.Val){
      pred.list[[v]]$v = SetValBias*pred.list[[v]]$v
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
  
  all_wts = array(data = NA, dim = c((maxit + 1), Unknown.ppp$n, n.sp))
  
  while(breakloop == 0)
  {
    niter = niter + 1
	
    #2 Compute predicted intensities
    
    pr_ppm_list = list()
    for (i in 1:nclust) {
      if(is.null(cov.bias))
      {
        pr_ppm_list[[i]] = predict(ppm_list[[i]], locations = Unknown.ppp)
      }
      else
      {
        pr_ppm_list[[i]] = predict(ppm_list[[i]], covariates = pred.list, locations = Unknown.ppp)
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
        addtosp.list[[i]] = (1:Unknown.ppp$n)
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
      if(num.add > Unknown.ppp$n/n.sp)
      {
        print("Impossible to add so many points, the highest possible number will be used instead")
        num.add = floor(Unknown.ppp$n/n.sp)
      }
      else
      {
        num.add = num.add
      }
      
      add_max = apply(test_wts, 2, sort, decreasing = TRUE)[num.add,]
      addtosp.list =list()
      for (i in 1:nclust) {
        addtosp.list[[i]] = if(anyNA(Unknown.ppp$x[test_wts[,i] >= add_max[i]]) == TRUE) integer() else which(test_wts[,i] >= add_max[i])
        i=i+1
      }
      
    }
    
    # lists and vectors needed in the next steps
    sp_aug.list = sp_wts.list = sp_aug_ppp.list = Q_aug.list = ppm.L.list = list()
    ppm.L.pred = ppm.pred.list = list()
    Dloglik = counts.sp = rep(NA, nclust)
    
    for (i in 1:nclust) {
      sp_aug.list[[i]] = data.frame(X = c(ppp_list[[i]]$x, Unknown.ppp$x[addtosp.list[[i]]]), Y = c(ppp_list[[i]]$y, Unknown.ppp$y[addtosp.list[[i]]])) # add unknown points to known points of species 1
      quad.xy = data.frame(quads.$x, quads.$y)
      names(quad.xy) = c("X", "Y")
      
      #5 update quadrature weights
      
      #scores for species with known label (weight =1) and for the new obs (test_wts)
      sp_wts.list[[i]] = scoreweights(sp_aug.list[[i]], quad.xy, scores = c(rep(1, ppp_list[[i]]$n), test_wts[addtosp.list[[i]], i])) # generate quad weights for augmented species 1
      
      # Augmented point patterns
      #win.   = win. # added because HPC wouldn't work
      
      sp_aug_ppp.list[[i]] = ppp(x =sp_aug.list[[i]]$X, y = sp_aug.list[[i]]$Y,
                                 window =winMix)
      
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
      
      if(is.null(rAreaInter))
      {
        ppm.L.list[[i]] = ppm(Q_aug.list[[i]], trend = ppmform, covariates = cov.list, 
                              gcontrol = list(epsilon = 1e-6, maxit = 100))
      }
      else
      {
        ppm.L.list[[i]] = ppm(Q_aug.list[[i]], trend = ppmform, covariates = cov.list,
                              AreaInter(rAreaInter),
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
      if (num.add > Unknown.ppp$n/n.sp)
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
      pred.loc[[l]] = predict(ppm_list[[l]], locations = quads.)
    }else{
      pred.loc[[l]] = predict(ppm_list[[l]], covariates = pred.list, 
                                        locations = quads.)
    }
  }
  
  fitaft.pred=list()
  # For predictions
  for (i in 1:nclust) {
  if(is.null(cov.bias))
  {
    fitaft.pred[[i]]  = predict(ppm_list[[i]], covariates=cov.list)
  }
  else
  {
    fitaft.pred[[i]]  = predict(ppm_list[[i]], covariates = pred.list)
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
              Weight.df = test_wts,
              fit.final = ppm_list,
              niter = niter, 
              ppm.pred.list = pr_ppm_list,
              ppm.pred.list.unk = pr_ppm_list.unk,
              sp_aug.list = sp_aug.list,
              sp_aug_ppp.list = sp_aug_ppp.list,
              delta_max = delta_max,
              delta_min = delta_min,
              delta_step = delta_step,
              classif = "Loop",
              pred.loc = pred.loc,
              fitaft.pred=fitaft.pred,
              sp_aug.list = sp_aug.list,
              Known.ppp = Known.ppp,
              Unknown.ppp =Unknown.ppp
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


IMSE = function(mu1, mu2, fun = "Else", mu.min = 1.e-5, rescale = TRUE)
{
  mu1.use = mu1
  mu2.use = mu2
  
  if (rescale == TRUE)
  {
    mu2.use = mu2.use*mean(mu1.use)/mean(mu2.use)
  }
  
  #mu1.use = mu1
  mu1.use[mu1.use < mu.min] = mu.min
  #mu2.use = mu2
  mu2.use[mu2.use < mu.min] = mu.min
  
  if (fun == "Else")
  {
    mu1.use = mu1.use
    mu2.use = mu2.use
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

corint = function(mu1, mu2, fun = "Else", method=c("pearson", "kendall", "spearman"), mu.min = 1.e-5)
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
  
  if (fun == "Else")
  {
    mu1.use = mu1.use
    mu2.use = mu2.use
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

RSS = function(weightmatrix, Unknown_labels)
{
  mark.cols = match(Unknown_labels, colnames(weightmatrix))
  correctweights = weightmatrix[cbind(seq_along(mark.cols), mark.cols)]
  RSSscore = sum((correctweights - 1)^2)
  RSSscore
}

###----------------- Accuracy
Accuracy = function(Known.ppp, New_weights, Unknown_labels){
  
  n.sp = length(unique(marks(Known.ppp)))
  
  W.max = apply(New_weights[(1:nrow(New_weights)),], 1, max)
  C.id = apply(New_weights[(1:nrow(New_weights)),], 1, which.max)
  
  indiv_testlab = as.data.frame(cbind(Unknown_labels, C.id, W.max))
  
  levels(indiv_testlab$Unknown_labels) <- c(unique(Unknown_labels))
  
  # New method for accuracy
  levels(indiv_testlab$Unknown_labels)
  levels(factor(indiv_testlab$C.id))
  for (i in 1:n.sp) {
    levels(indiv_testlab$Unknown_labels)[levels(indiv_testlab$Unknown_labels)==paste("Sp", i, sep = "")] <- paste(i, sep = "")
    
  }
  
  acc_tab <- table(indiv_testlab$Unknown_labels,
                   factor(indiv_testlab$C.id, levels=1:n.sp))
  
  #dim(acc_tab)
  
  #CM.acc = confusionMatrix(acc_tab)
  m.acc= acc_tab
  
  # some usuful calc
  n = sum(m.acc) # number of instances
  nc = nrow(m.acc) # number of classes
  diag = diag(m.acc) # number of correctly classified instances per class 
  rowsums = apply(m.acc, 1, sum) # number of instances per class
  colsums = apply(m.acc, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  
  # accuracy measure
  accuracy = sum(diag) / n 
  accuracy 
  
}


Perffunc = function(fit, sp.int, fun = "Else", rescale = TRUE, 
                    method=c("pearson", "kendall", "spearman"), 
                    mu.min = 1.e-5, Known.ppp., Unknown_labels., pf = c(NULL)){
  
  # Set up otpion for initial weights  
  n.sp = length(unique(marks(Known.ppp)))
  
  datappp = superimpose.ppp(Known.ppp, Unknown.ppp)
  pred = fit$pred.loc
  New_weights = fit$Newpts_w
  
  if(fit$classif=="Loop"){
    colnames(New_weights)= c(unique(Unknown_labels))
  }
  
  ###----------------- IMSE
  
  IMSE = function(mu1, mu2, fun = "Else", mu.min = 1.e-5, rescale = TRUE)
  {
    mu1.use = mu1
    mu2.use = mu2
    
    if (rescale == TRUE)
    {
      mu2.use = mu2.use*mean(mu1.use)/mean(mu2.use)
    }
    
    #mu1.use = mu1
    mu1.use[mu1.use < mu.min] = mu.min
    #mu2.use = mu2
    mu2.use[mu2.use < mu.min] = mu.min
    
    if (fun == "Else")
    {
      mu1.use = mu1.use
      mu2.use = mu2.use
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
  
  corint = function(mu1, mu2, fun = "Else", method=c("pearson", "kendall", "spearman"), mu.min = 1.e-5)
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
    
    if (fun == "Else")
    {
      mu1.use = mu1.use
      mu2.use = mu2.use
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
  
  RSS = function(weightmatrix, Unknown_labels)
  {
    mark.cols = match(Unknown_labels, colnames(weightmatrix))
    correctweights = weightmatrix[cbind(seq_along(mark.cols), mark.cols)]
    RSSscore = sum((correctweights - 1)^2)
    RSSscore
  }
  
  ###----------------- Accuracy
  
  Accuracy = function(Known.ppp, New_weights, Unknown_labels){
    
    n.sp = length(unique(marks(Known.ppp)))
    
    W.max = apply(New_weights[(1:nrow(New_weights)),], 1, max)
    C.id = apply(New_weights[(1:nrow(New_weights)),], 1, which.max)
    
    indiv_testlab = as.data.frame(cbind(Unknown_labels, C.id, W.max))
    
    levels(indiv_testlab$Unknown_labels) <- c(unique(Unknown_labels))
    
    # New method for accuracy
    levels(indiv_testlab$Unknown_labels)
    levels(as.factor(indiv_testlab$C.id))
    
      for (i in 1:n.sp) {
        indiv_testlab$Unknown_labels[levels(indiv_testlab$Unknown_labels)==paste("Sp", i, sep = "")] <- paste(i, sep = "")
        
      }

    
    
    acc_tab <- table(factor(indiv_testlab$Unknown_labels, levels=1:n.sp), 
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
    sp.predlist = IMSEs = list()
    
    if (fit$LoopM == TRUE){
      for (l in 1:n.sp) {
        sp.predlist[[l]] = as.vector(t(pred[[l]]))
        IMSEs[[l]] = IMSE(sp.int[[l]], sp.predlist[[l]])
      }
    }else{
        for (l in 1:n.sp) {
          sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
          IMSEs[[l]] = IMSE(sp.int[[l]], sp.predlist[[l]])
        }
    }
    
    IMSEscore = sum(unlist(IMSEs))
    
  }
  
  if(any(pf)=="corint"){
    sp.predlist = sumcor1l = sumcor2l = list()
    for (l in 1:n.sp) {
      sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
      sumcor1l = corint(sp.int[[l]], sp.predlist[[l]], method="pearson")
      sumcor2l = corint(sp.int[[l]], sp.predlist[[l]], method="kendall")
    }
    sumcor1 = sum(unlist(sumcor1l))
    sumcor2 = sum(unlist(sumcor2l))
  }
  
  if(any(pf)=="mRSS"){
    #RSS1 = RSS(New_weights, Unknown_labels)
    meanRSS = RSS(New_weights, Unknown_labels)/length(Unknown_labels)
  }
  
  if(any(pf)=="Acc"){
    accmat = Accuracy(Known.ppp, New_weights, Unknown_labels)
  }
  
  if(is.null(pf)==TRUE){
    accmat = Accuracy(Known.ppp, New_weights, Unknown_labels)
    #RSSscore = RSS(New_weights, Unknown_labels)
    meanRSS = RSS(New_weights, Unknown_labels)/length(Unknown_labels)
    
    sp.predlist = IMSEs = sumcor1l = sumcor2l = list()
    
    if (fit$classif=="Loop"){
      for (l in 1:n.sp) {
        sp.predlist[[l]] = as.vector(t(pred[[l]]))
        IMSEs[[l]] = IMSE(sp.int[[l]], sp.predlist[[l]])
      }
    }else{
        for (l in 1:n.sp) {
          sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
          IMSEs[[l]] = IMSE(sp.int[[l]], sp.predlist[[l]])
        }
    } 
    IMSEscore = sum(unlist(IMSEs))
    
    for (l in 1:n.sp) {
      #sp.predlist[[l]] = as.vector(t(pred[[l]]$v))
      sumcor1l[[l]] = corint(sp.int[[l]], sp.predlist[[l]], method="pearson")
      sumcor2l[[l]] = corint(sp.int[[l]], sp.predlist[[l]], method="kendall")
    }
    sumcor1 = sum(unlist(sumcor1l))
    sumcor2 = sum(unlist(sumcor2l))
    
    
  }
  
  return(list(accmat=accmat, meanRSS=meanRSS, sumIMSE=IMSEscore, sumcor1=sumcor1, sumcor2=sumcor2))
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

Testsims = function(hidepct, n.sims, sp_sim.list, k = k, ks=ks, 
                    nstart=nstart, quadsenv, cov.bias=NULL, SetValBias=NULL,
                    rAreaInter=NULL, delta_max=delta_max, delta_min=delta_min, 
                    delta_step =delta_step, num.add = num.add)
{
  # prepare quadrature ppp
  winMix   = owin(xrange = c(min(quadsenv$X)-0.5, max(quadsenv$X)+0.5), 
                  yrange = c(min(quadsenv$Y)-0.5, max(quadsenv$Y)+0.5))
  quads. = ppp(quadsenv$X, quadsenv$Y, window = winMix)
  Qmask = makeMask(quads.)
  
  
  RSSknn = meanRSSknn = IMSEknn = IMSEknn2 = IMSEknn3 = RSSkmeans = meanRSSkmeans = IMSEkmeans =  IMSEkmeans2 = IMSEkmeans3 = 
    RSSkps = meanRSSkps = IMSEkps = IMSEkps2 = IMSEkps3 = 
    RSSrand = meanRSSrand = IMSErand = IMSErand2 = IMSErand3 = RSSequal = meanRSSequal = IMSEequal =  IMSEequal2 = 
    IMSEequal3 = RSSCF = meanRSSCF = IMSECF = IMSECF2 = IMSECF3 = 
    RSSindiv = meanRSSindiv = IMSEindiv = IMSEindiv2 = IMSEindiv3 = RSSLoopT = meanRSSLoopT = IMSELoopT = IMSELoopT2 = IMSELoopT3 =
    RSSLoopE = meanRSSLoopE = IMSELoopE = IMSELoopE2 = IMSELoopE3 = RSSLoopA = meanRSSLoopA = IMSELoopA = IMSELoopA2 = IMSELoopA3 =
    sumcorknn1 = sumcorkmeans1 = sumcorrand1 = sumcorequal1 = sumcorindiv1 = sumcorLoopE1 =
    sumcorLoopA1 = sumcorLoopT1 = sumcorkps1 = sumcorCF1 = 
    sumcorknn2 = sumcorkmeans2 = sumcorrand2 = sumcorequal2 = sumcorindiv2 = sumcorLoopE2 =
    sumcorLoopA2 = sumcorLoopT2 = sumcorknn3 = sumcorkmeans3 = sumcorrand3 = sumcorequal3 = sumcorindiv3 = sumcorLoopE3 =
    sumcorLoopA3 = sumcorLoopT3 = sumcorkps2 = sumcorCF2 = sumcorkps3 = sumcorCF3 = matrix(NA, n.sims, length(hidepct))
  accmatknn = accmatkmeans = accmatrand = accmatequal = accmatindiv = accmatLoop = 
    accmatLoopA = accmatLoopT = accmatLoopE = accmatkps = accmatCF = matrix(NA, n.sims, length(hidepct))
  
  
  knnpred = kmeanspred = randpred = equalpred = LoopApred = LoopTpred = LoopEpred =
    indivpred = kpspred = CFpred = array(NA, c(quads.$n, 3, n.sims, length(hidepct)))
  
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
      
      Unknown.pp = ppp(x = c(unlist(coordtestx.list)), 
                     y = c(unlist(coordtesty.list)), window = win,
                     marks = c(unlist(markshide.list)))
      
      Unknown.ppp = ppp(x = c(unlist(coordtestx.list)), 
                      y = c(unlist(coordtesty.list)), window = win,
                      marks = c(rep("Unknown", Unknown.ppp$n)))
      
      Unknown_labels = as.vector(unlist(markstest.list))
      
      Known.ppp = ppp(x = c(unlist(coordsubx.list)), 
                     y = c(unlist(coordsuby.list)), window = win,
                     marks = c(unlist(marksub.list)))
      
      
      datappp = superimpose.ppp(Known.ppp, Unknown.ppp)
      
      n.sp = length(unique(marks(Known.ppp)))
      
      
      ###
      #  Mixture model
      ###---
      simknn = ppmMixEngine(Known.ppp, Unknown.ppp, quadsenv = Quadmat,
                            initweights = "knn", #n.sp=n.sp, 
                            k=k, ks=ks, nstart=nstart, ppmform = ppmform, 
                            cov.bias = cov.bias, SetValBias = SetValBias, rAreaInter = rAreaInter,
                            verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simkmeans = ppmMixEngine(Known.ppp, Unknown.ppp, quadsenv = Quadmat,
                               initweights = "kmeans", #n.sp=n.sp,
                               k=k, ks=ks, nstart=nstart, ppmform = ppmform, 
                               cov.bias = cov.bias, SetValBias = SetValBias, rAreaInter = rAreaInter,
                               verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simrandom = ppmMixEngine(Known.ppp, Unknown.ppp, quadsenv = Quadmat,
                               initweights = "random", #n.sp=n.sp,
                               k=k, ks=ks, nstart=NULL, ppmform = ppmform, 
                               cov.bias = cov.bias, SetValBias = SetValBias, rAreaInter = rAreaInter,
                               verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simCF = ppmMixEngine(Known.ppp, Unknown.ppp, quadsenv = Quadmat,
                           initweights = "CoinF", #n.sp=n.sp,
                           k=NULL, ks=NULL, nstart=NULL, ppmform = ppmform,
                           cov.bias = cov.bias, SetValBias = SetValBias, rAreaInter = rAreaInter,
                           verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      simkps = ppmMixEngine(Known.ppp, Unknown.ppp, quadsenv = Quadmat, 
                            initweights = "kps", #n.sp=n.sp,
                            k=NULL, ks=ks, nstart=nstart, ppmform = ppmform,
                            cov.bias = cov.bias, SetValBias = SetValBias, rAreaInter = rAreaInter,
                            verbose = TRUE, tol = 0.000001, maxit = 50, plots = FALSE)
      
      # for performance measures
      knn_weights = simknn$New_weights[((Known.ppp$n)+1):nrow(simknn$New_weights),]
      pred.knn    = knn_weights
      W.knn[[i]][[j]] = knn_weights
      
      kmeans_weights = simkmeans$New_weights[((Known.ppp$n)+1):nrow(simkmeans$New_weights),]
      pred.kmeans    = kmeans_weights
      W.kmeans[[i]][[j]] = kmeans_weights
      
      random_weights = simrandom$New_weights[((Known.ppp$n)+1):nrow(simrandom$New_weights),]
      pred.random    = random_weights
      W.rand[[i]][[j]] = random_weights
      
      kps_weights = simkps$New_weights[((Known.ppp$n)+1):nrow(simkps$New_weights),]
      pred.kps    = kps_weights
      W.kps[[i]][[j]] = kps_weights
      
      #equal_weights = simequal$New_weights[((Known.ppp$n)+1):nrow(simequal$New_weights),]
      #pred.equal    = equal_weights
      
      CF_weights = simCF$New_weights[((Known.ppp$n)+1):nrow(simCF$New_weights),]
      pred.CF    = CF_weights
      W.CF[[i]][[j]] = CF_weights
      
      
      accmatknn[j, i] = Accuracy(Known.ppp, knn_weights, Unknown_labels)
      RSSknn[j, i] = RSS(pred.knn, Unknown_labels)
      meanRSSknn[j, i] = RSS(pred.knn, Unknown_labels)/length(Unknown_labels)
      
      accmatkmeans[j, i] = Accuracy(Known.ppp, kmeans_weights, Unknown_labels)
      RSSkmeans[j, i] = RSS(pred.kmeans, Unknown_labels)
      meanRSSkmeans[j, i] = RSS(pred.kmeans, Unknown_labels)/length(Unknown_labels)
      
      accmatrand[j, i] = Accuracy(Known.ppp, random_weights, Unknown_labels)
      RSSrand[j, i] = RSS(pred.random, Unknown_labels)
      meanRSSrand[j, i] = RSS(pred.random, Unknown_labels)/length(Unknown_labels)
      
      #accmatequal[j, i] = Accuracy(Known.ppp, equal_weights, Unknown_labels, n.sp)
      #RSSequal[j, i] = RSS(pred.equal, Unknown_labels)
      #meanRSSequal[j, i] = RSS(pred.equal, Unknown_labels)/length(Unknown_labels)
      
      accmatCF[j, i] = Accuracy(Known.ppp, CF_weights, Unknown_labels)
      RSSCF[j, i] = RSS(pred.CF, Unknown_labels)
      meanRSSCF[j, i] = RSS(pred.CF, Unknown_labels)/length(Unknown_labels)
      
      accmatkps[j, i] = Accuracy(Known.ppp, kps_weights, Unknown_labels)
      RSSkps[j, i] = RSS(pred.kps, Unknown_labels)
      meanRSSkps[j, i] = RSS(pred.kps, Unknown_labels)/length(Unknown_labels)
      
      #--
      if(is.null(cov.bias)){
        pred.knn = predict(simknn$fit.final, locations = quads.)
      }else{
        pred.knn = predict(simknn$fit.final, covariates = pred.list, locations = quads.)
      }
      
      sp.predlist.knn = IMSEknn.list = IMSEknn2.list = IMSEknn3.list = sumcorknn1.list = sumcorknn2.list = sumcorknn3.list = list()
      for (l in 1:n.sp) {
        sp.predlist.knn[[l]] = as.vector(t(pred.knn[[l]]))
        
        IMSEknn.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.knn[[l]], fun="log")
        IMSEknn2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.knn[[l]])
        IMSEknn3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.knn[[l]], fun="log", mu.min=1e-10)
        sumcorknn1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.knn[[l]], method="pearson", fun="log")
        sumcorknn2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.knn[[l]], method="pearson")
        sumcorknn3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.knn[[l]], fun="log", method="pearson", mu.min=1e-10)
        
      }
      
      IMSEknn[j, i] = sum(unlist(IMSEknn.list))
      IMSEknn2[j, i] = sum(unlist(IMSEknn2.list))
      IMSEknn3[j, i] = sum(unlist(IMSEknn3.list))
      sumcorknn1[j, i] = sum(unlist(sumcorknn1.list))
      sumcorknn2[j, i] = sum(unlist(sumcorknn2.list))
      sumcorknn3[j, i] = sum(unlist(sumcorknn3.list))
      
      
      #--
      if(is.null(cov.bias)){
        pred.kmeans = predict(simkmeans$fit.final, locations = quads.)
      }else{
        pred.kmeans = predict(simkmeans$fit.final, covariates = pred.list, locations = quads.)
      }
      
      sp.predlist.kmeans = IMSEkmeans.list = IMSEkmeans2.list = IMSEkmeans3.list = sumcorkmeans1.list = sumcorkmeans2.list = sumcorkmeans3.list = list()
      for (l in 1:n.sp) {
        sp.predlist.kmeans[[l]] = as.vector(t(pred.kmeans[[l]]))
        
        IMSEkmeans.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.kmeans[[l]], fun="log")
        IMSEkmeans2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.kmeans[[l]])
        IMSEkmeans3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.kmeans[[l]], fun="log", mu.min=1e-10)
        sumcorkmeans1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.kmeans[[l]], method="pearson", fun="log")
        sumcorkmeans2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.kmeans[[l]], method="pearson")
        sumcorkmeans3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.kmeans[[l]], fun="log", method="pearson", mu.min=1e-10)
        
      }
      
      IMSEkmeans[j, i] = sum(unlist(IMSEkmeans.list))
      IMSEkmeans2[j, i] = sum(unlist(IMSEkmeans2.list))
      IMSEkmeans3[j, i] = sum(unlist(IMSEkmeans3.list))
      sumcorkmeans1[j, i] = sum(unlist(sumcorkmeans1.list))
      sumcorkmeans2[j, i] = sum(unlist(sumcorkmeans2.list))
      sumcorkmeans3[j, i] = sum(unlist(sumcorkmeans3.list))
      
      #--
      if(is.null(cov.bias)){
        pred.random  = predict(simrandom$fit.final, locations = quads.)
      }else{
        pred.random = predict(simrandom$fit.final, covariates = pred.list, locations = quads.)
      }
      
      sp.predlist.rand = IMSErand.list = IMSErand2.list = IMSErand3.list = sumcorrand1.list = sumcorrand2.list = sumcorrand3.list = list()
      for (l in 1:n.sp) {
        sp.predlist.rand[[l]] = as.vector(t(pred.random[[l]]))
        
        IMSErand.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.rand[[l]], fun="log")
        IMSErand2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.rand[[l]])
        IMSErand3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.rand[[l]], fun="log", mu.min=1e-10)
        sumcorrand1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.rand[[l]], method="pearson", fun="log")
        sumcorrand2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.rand[[l]], method="pearson")
        sumcorrand3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.rand[[l]], fun="log", method="pearson", mu.min=1e-10)
        
      }
      
      IMSErand[j, i] = sum(unlist(IMSErand.list))
      IMSErand2[j, i] = sum(unlist(IMSErand2.list))
      IMSErand3[j, i] = sum(unlist(IMSErand3.list))
      sumcorrand1[j, i] = sum(unlist(sumcorrand1.list))
      sumcorrand2[j, i] = sum(unlist(sumcorrand2.list))
      sumcorrand3[j, i] = sum(unlist(sumcorrand3.list))
      
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
        pred.CF  = predict(simCF$fit.final, locations = quads.)
      }else{
        pred.CF = predict(simCF$fit.final, covariates = pred.list, 
                          locations = quads.)
      }
      
      sp.predlist.CF = IMSECF.list = IMSECF2.list = IMSECF3.list = sumcorCF1.list = sumcorCF2.list = sumcorCF3.list = list()
      for (l in 1:n.sp) {
        sp.predlist.CF[[l]] = as.vector(t(pred.CF[[l]]))
        
        IMSECF.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.CF[[l]], fun="log")
        IMSECF2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.CF[[l]])
        IMSECF3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.CF[[l]], fun="log", mu.min=1e-10)
        sumcorCF1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.CF[[l]], method="pearson", fun="log")
        sumcorCF2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.CF[[l]], method="pearson")
        sumcorCF3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.CF[[l]], method="pearson", fun="log", mu.min=1e-10)
        
      }
      
      IMSECF[j, i] = sum(unlist(IMSECF.list))
      IMSECF2[j, i] = sum(unlist(IMSECF2.list))
      IMSECF3[j, i] = sum(unlist(IMSECF3.list))
      sumcorCF1[j, i] = sum(unlist(sumcorCF1.list))
      sumcorCF2[j, i] = sum(unlist(sumcorCF2.list))
      sumcorCF3[j, i] = sum(unlist(sumcorCF3.list))
      
      #--
      if(is.null(cov.bias)){
        pred.kps = predict(simkps$fit.final, locations = quads.)
      }else{
        pred.kps = predict(simkps$fit.final, covariates = pred.list, 
                           locations = quads.)
      }
      
      sp.predlist.kps = IMSEkps.list = IMSEkps2.list = IMSEkps3.list = sumcorkps1.list = sumcorkps2.list = sumcorkps3.list = list()
      for (l in 1:n.sp) {
        sp.predlist.kps[[l]] = as.vector(t(pred.kps[[l]]))
        
        IMSEkps.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.kps[[l]], fun="log")
        IMSEkps2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.kps[[l]])
        IMSEkps3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.kps[[l]], fun="log", mu.min=1e-10)
        sumcorkps1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.kps[[l]], method="pearson", fun="log")
        sumcorkps2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.kps[[l]], method="pearson")
        sumcorkps3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.kps[[l]], method="pearson", fun="log", mu.min=1e-10)
        
      }
      
      IMSEkps[j, i] = sum(unlist(IMSEkps.list))
      IMSEkps2[j, i] = sum(unlist(IMSEkps2.list))
      IMSEkps3[j, i] = sum(unlist(IMSEkps3.list))
      sumcorkps1[j, i] = sum(unlist(sumcorkps1.list))
      sumcorkps2[j, i] = sum(unlist(sumcorkps2.list))
      sumcorkps3[j, i] = sum(unlist(sumcorkps3.list))
      
      
      # for intensity plots
      for (l in 1:n.sp) {
        knnpred[,l,j,i] = as.matrix(sp.predlist.knn[[l]])
        kmeanspred[,l,j,i] = as.matrix(sp.predlist.kmeans[[l]])
        randpred[,l,j,i] = as.matrix(sp.predlist.rand[[l]])
        #equalpred[,,j,i] = as.matrix(unlist(sp.predlist.equal))
        CFpred[,l,j,i] = as.matrix(sp.predlist.CF[[l]])
        kpspred[,l,j,i] = as.matrix(sp.predlist.kps[[l]])
      }
      
      
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
      
      simLoopA = ppmLoopEngine(Known.ppp, Unknown.ppp, addpt = "LoopA", quadsenv=quadsenv,
                               ppmform=ppmform, delta_max=NULL, delta_min=NULL, delta_step=NULL, num.add = NULL,
                               cov.bias=NULL, SetValBias =NULL, rAreaInter=NULL, maxit = 50,
                               tol=0.000001, verbose = TRUE, plots = FALSE)
      
      simLoopT = ppmLoopEngine(Known.ppp, Unknown.ppp, addpt = "LoopT", quadsenv=quadsenv,
                               ppmform= ppmform, delta_max=delta_max, delta_min=delta_min, delta_step=delta_step, num.add = NULL,
                               cov.bias=NULL, SetValBias =NULL, rAreaInter=NULL, maxit = 50,
                               tol=0.000001, verbose = TRUE, plots = FALSE)
      
      simLoopE = ppmLoopEngine(Known.ppp, Unknown.ppp, addpt = "LoopE", quadsenv=quadsenv,
                               ppmform= ppmform, delta_max=NULL, delta_min=NULL, delta_step=NULL, num.add = num.add,
                               cov.bias=NULL, SetValBias =NULL, rAreaInter=NULL, maxit = 50,
                               tol=0.000001, verbose = TRUE, plots = FALSE)
      
      
      # for performance measures
      
      LoopA_weights = simLoopA$New_weights
      pred.LoopA    = LoopA_weights
      colnames(pred.LoopA)= c(unique(Unknown_labels))
      W.LA[[i]][[j]] = LoopA_weights
      
      LoopT_weights = simLoopT$New_weights
      pred.LoopT    = LoopT_weights
      colnames(pred.LoopT)= c(unique(Unknown_labels))
      W.LT[[i]][[j]] = LoopT_weights
      
      LoopE_weights = simLoopE$New_weights
      pred.LoopE    = LoopE_weights
      colnames(pred.LoopE)= c(unique(Unknown_labels))
      W.LE[[i]][[j]] = LoopE_weights
      
      #
      accmatLoopA[j, i] = Accuracy(Known.ppp, LoopA_weights, Unknown_labels)
      RSSLoopA[j, i] = RSS(pred.LoopA, Unknown_labels)
      meanRSSLoopA[j, i] = RSS(pred.LoopA, Unknown_labels)/length(Unknown_labels)
      
      accmatLoopT[j, i] = Accuracy(Known.ppp, LoopT_weights, Unknown_labels)
      RSSLoopT[j, i] = RSS(pred.LoopT, Unknown_labels)
      meanRSSLoopT[j, i] = RSS(pred.LoopT, Unknown_labels)/length(Unknown_labels)
      
      accmatLoopE[j, i] = Accuracy(Known.ppp, LoopE_weights, Unknown_labels)
      RSSLoopE[j, i] = RSS(pred.LoopE, Unknown_labels)
      meanRSSLoopE[j, i] = RSS(pred.LoopE, Unknown_labels)/length(Unknown_labels)
      
      #--
      pr_quad_ppmlist.LA = pr_quad_ppmlist.LT = pr_quad_ppmlist.LE = list()
      for (l in 1:n.sp) {
        if(is.null(cov.bias)){
          pr_quad_ppmlist.LA[[l]] = predict(simLoopA$ppm_list[[l]], locations = quads.)
        }else{
          pr_quad_ppmlist.LA[[l]] = predict(simLoopA$ppm_list[[l]], covariates = pred.list, locations = quads.)
        }
        
        if(is.null(cov.bias)){
          pr_quad_ppmlist.LT[[l]] = predict(simLoopT$ppm_list[[l]], locations = quads.)
        }else{
          pr_quad_ppmlist.LT[[l]] = predict(simLoopT$ppm_list[[l]], covariates = pred.list, locations = quads.)
        }
        
        if(is.null(cov.bias)){
          pr_quad_ppmlist.LE[[l]] = predict(simLoopE$ppm_list[[l]], locations = quads.)
        }else{
          pr_quad_ppmlist.LE[[l]] = predict(simLoopE$ppm_list[[l]], covariates = pred.list, locations = quads.)
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
      
      sp.predlist.LA = IMSELA.list = IMSELA2.list = IMSELA3.list = sumcorLA1.list = sumcorLA2.list = sumcorLA3.list = list()
      sp.predlist.LT = IMSELT.list = IMSELT2.list = IMSELT3.list = sumcorLT1.list = sumcorLT2.list = sumcorLT3.list = list()
      sp.predlist.LE = IMSELE.list = IMSELE2.list = IMSELE3.list = sumcorLE1.list = sumcorLE2.list = sumcorLE3.list = list()
      for (l in 1:n.sp) {
        sp.predlist.LA[[l]] = as.vector(t(pr_quad_ppmlist.LA[[l]]))
        
        IMSELA.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LA[[l]], fun="log")
        IMSELA2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LA[[l]])
        IMSELA3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LA[[l]], fun="log", mu.min=1e-10)
        sumcorLA1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LA[[l]], method="pearson", fun="log")
        sumcorLA2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LA[[l]], method="pearson")
        sumcorLA3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LA[[l]], method="pearson", fun="log", mu.min=1e-10)
        
        sp.predlist.LT[[l]] = as.vector(t(pr_quad_ppmlist.LT[[l]]))
        
        IMSELT.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LT[[l]], fun="log")
        IMSELT2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LT[[l]])
        IMSELT3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LT[[l]], fun="log", mu.min=1e-10)
        sumcorLT1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LT[[l]], method="pearson", fun="log")
        sumcorLT2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LT[[l]], method="pearson")
        sumcorLT3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LT[[l]], method="pearson", fun="log", mu.min=1e-10)
        
        sp.predlist.LE[[l]] = as.vector(t(pr_quad_ppmlist.LE[[l]]))
        
        IMSELE.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LE[[l]], fun="log")
        IMSELE2.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LE[[l]])
        IMSELE3.list[[l]] = IMSE(sp_int.list[[l]], sp.predlist.LE[[l]], fun="log", mu.min=1e-10)
        sumcorLE1.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LE[[l]], method="pearson", fun="log")
        sumcorLE2.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LE[[l]], method="pearson")
        sumcorLE3.list[[l]] = corint(sp_int.list[[l]], sp.predlist.LE[[l]], method="pearson", fun="log", mu.min=1e-10)
      }
      
      IMSELoopA[j, i] = sum(unlist(IMSELA.list))
      IMSELoopA2[j, i] = sum(unlist(IMSELA2.list))
      IMSELoopA3[j, i] = sum(unlist(IMSELA3.list))
      sumcorLoopA1[j, i] = sum(unlist(sumcorLA1.list))
      sumcorLoopA2[j, i] = sum(unlist(sumcorLA2.list))
      sumcorLoopA3[j, i] = sum(unlist(sumcorLA3.list))
      
      IMSELoopT[j, i] = sum(unlist(IMSELT.list))
      IMSELoopT2[j, i] = sum(unlist(IMSELT2.list))
      IMSELoopT3[j, i] = sum(unlist(IMSELT3.list))
      sumcorLoopT1[j, i] = sum(unlist(sumcorLT1.list))
      sumcorLoopT2[j, i] = sum(unlist(sumcorLT2.list))
      sumcorLoopT3[j, i] = sum(unlist(sumcorLT3.list))
      
      IMSELoopE[j, i] = sum(unlist(IMSELE.list))
      IMSELoopE2[j, i] = sum(unlist(IMSELE2.list))
      IMSELoopE3[j, i] = sum(unlist(IMSELE3.list))
      sumcorLoopE1[j, i] = sum(unlist(sumcorLE1.list))
      sumcorLoopE2[j, i] = sum(unlist(sumcorLE2.list))
      sumcorLoopE3[j, i] = sum(unlist(sumcorLE3.list))
      
      
      
      ####
      # Individual PPMs
      ###---
      
      datamarks = marks(datappp)
      uniquemarks = unique(datamarks)
      colmarks = uniquemarks
      unknown = datamarks == "Unknown"
      names.mark  = colmarks[-which(colmarks == "Unknown")]
       
      
      Qind = ppm_spind = pr_sp_ppmind = list()
      for (l in 1:n.sp) {
        #specie separetely
        Qind[[l]] = quadscheme(data = sp_sub.list[[l]], dummy = quads., method = "grid", ntile = c(101, 101), npix = c(101, 101))
        
        # Fit Poisson PPMs
        if(is.null(rAreaInter)){
          ppm_spind[[l]]  = ppm(Qind[[l]], trend = ppmform, covariates = cov.list)
          
        }else{
          ppm_spind[[l]]  = ppm(Qind[[l]], trend = ppmform, covariates = cov.list, AreaInter(rAreaInter))
        }
        
        # Predict intensity at locations where labels are hidden
        ## former label from sp1 tested for the 3PPMs
        if(is.null(cov.bias)){
          pr_sp_ppmind[[l]] = predict(ppm_spind[[l]], locations = Unknown.ppp)
        }else{
          pr_sp_ppmind[[l]] = predict(ppm_spind[[l]], covariates = pred.list, locations = Unknown.ppp)
        }
      }
      
      
      sp_predsind = data.frame(matrix(unlist(pr_sp_ppmind),
                                      nrow=length(pr_sp_ppmind[[1]]), byrow=F))
      ind_wts  = sp_predsind/apply(sp_predsind, 1, sum)
      max_predind = apply(sp_predsind, 1, which.max)
      max_predind.vec = as.vector(table(max_predind))
      
      # calculate accuracy in the same way that we did for the mixture models
      checkind.id = as.vector(max_predind)
      checkind.data = as.data.frame(cbind(Unknown_labels, checkind.id))
      
      #to deal with some labels not be present in some iterations
      
      testindlab = as.data.frame(checkind.data$Unknown_labels)
      checkindlab = as.data.frame(checkind.data$check.id)
      
      checkindD = as.data.frame(checkind.data)
      
      m.acc.indiv = table(checkindD$Unknown_labels, checkindD$checkind.id)
      
      
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
      
      RSSindiv[j, i] = RSS(testind_wts, Unknown_labels)
      meanRSSindiv[j, i] = RSS(testind_wts, Unknown_labels)/length(Unknown_labels)
      
      pr_quad_ppmindlist = list()
      for (l in 1:n.sp) {
        if(is.null(cov.bias)){
          pr_quad_ppmindlist[[l]] = predict(ppm_spind[[l]], locations = quads.)
        }else{
          pr_quad_ppmindlist[[l]] = predict(ppm_spind[[l]], covariates = pred.list, locations = quads.)
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
      
      sp.predindlist = IMSEindiv.list = IMSEindiv3.list = IMSEindiv2.list = sumcorindiv1.list = sumcorindiv2.list = sumcorindiv3.list = list()
      for (l in 1:n.sp) {
        sp.predindlist[[l]] = as.vector(t(pr_quad_ppmindlist[[l]]))
        IMSEindiv.list[[l]] = IMSE(sp_int.list[[l]], sp.predindlist[[l]], fun="log")
        IMSEindiv2.list[[l]] = IMSE(sp_int.list[[l]], sp.predindlist[[l]])
        IMSEindiv3.list[[l]] = IMSE(sp_int.list[[l]], sp.predindlist[[l]], fun="log", mu.min=1e-10)
        sumcorindiv1.list[[l]] = corint(sp_int.list[[l]], sp.predindlist[[l]], method="pearson", fun="log")
        sumcorindiv2.list[[l]] = corint(sp_int.list[[l]], sp.predindlist[[l]], method="pearson")
        sumcorindiv3.list[[l]] = corint(sp_int.list[[l]], sp.predindlist[[l]], method="pearson", fun="log", mu.min=1e-10)
        
      }
      
      IMSEindiv[j, i] = sum(unlist(IMSEindiv.list))
      IMSEindiv2[j, i] = sum(unlist(IMSEindiv2.list))
      IMSEindiv3[j, i] = sum(unlist(IMSEindiv3.list))
      sumcorindiv1[j, i] = sum(unlist(sumcorindiv1.list))
      sumcorindiv2[j, i] = sum(unlist(sumcorindiv2.list))
      sumcorindiv3[j, i] = sum(unlist(sumcorindiv3.list))
      
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
      
      if(length(hidepct)==1){
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
  return(list(RSSknn = RSSknn, meanRSSknn = meanRSSknn, IMSEknn = IMSEknn, IMSEknn2 = IMSEknn2, IMSEknn3 = IMSEknn3, sumcorknn1 = sumcorknn1, sumcorknn2 = sumcorknn2, sumcorknn3 = sumcorknn3,
              RSSkmeans = RSSkmeans, meanRSSkmeans = meanRSSkmeans, IMSEkmeans = IMSEkmeans, IMSEkmeans2 = IMSEkmeans2, IMSEkmeans3 = IMSEkmeans3, sumcorkmeans1 = sumcorkmeans1, sumcorkmeans2 = sumcorkmeans2, sumcorkmeans3 = sumcorkmeans3, 
              RSSrand = RSSrand, meanRSSrand = meanRSSrand, IMSErand = IMSErand, IMSErand2 = IMSErand2, IMSErand3 = IMSErand3, sumcorrand1 = sumcorrand1, sumcorrand2 = sumcorrand2, sumcorrand3 = sumcorrand3,
              #RSSequal = RSSequal, meanRSSequal = meanRSSequal, IMSEequal = IMSEequal, sumcorequal2 = sumcorequal2, sumcorequal2 = sumcorequal2, 
              RSSkps = RSSkps, meanRSSkps = meanRSSkps, IMSEkps = IMSEkps, IMSEkps2 = IMSEkps2, IMSEkps3 = IMSEkps3, sumcorkps1 = sumcorkps1, sumcorkps2 = sumcorkps2, sumcorkps3 = sumcorkps3,
              RSSCF = RSSCF, meanRSSCF = meanRSSCF, IMSECF = IMSECF, IMSECF2 = IMSECF2, IMSECF3 = IMSECF3, sumcorCF1 = sumcorCF1, sumcorCF2 = sumcorCF2, sumcorCF3 = sumcorCF3,
              RSSindiv = RSSindiv, meanRSSindiv = meanRSSindiv, IMSEindiv = IMSEindiv, IMSEindiv2 = IMSEindiv2, IMSEindiv3 = IMSEindiv3, sumcorindiv1 = sumcorindiv1, sumcorindiv2 = sumcorindiv2, sumcorindiv3 = sumcorindiv3, 
              RSSLoopA = RSSLoopA, meanRSSLoopA = meanRSSLoopA, IMSELoopA = IMSELoopA, IMSELoopA2 = IMSELoopA2, IMSELoopA3 = IMSELoopA3, sumcorLoopA1 = sumcorLoopA1, sumcorLoopA2 = sumcorLoopA2, sumcorLoopA3 = sumcorLoopA3,
              RSSLoopT = RSSLoopT, meanRSSLoopT = meanRSSLoopT, IMSELoopT = IMSELoopT, IMSELoopT2 = IMSELoopT2, IMSELoopT3 = IMSELoopT3, sumcorLoopT1 = sumcorLoopT1, sumcorLoopT2 = sumcorLoopT2, sumcorLoopT3 = sumcorLoopT3, 
              RSSLoopE = RSSLoopE, meanRSSLoopE = meanRSSLoopE, IMSELoopE = IMSELoopE, IMSELoopE2 = IMSELoopE2, IMSELoopE3 = IMSELoopE3, sumcorLoopE1 = sumcorLoopE1, sumcorLoopE2 = sumcorLoopE2, sumcorLoopE3 = sumcorLoopE3, 
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

createConfusionMatrix <- function(act, pred) {
  pred <- pred[order(act)]
  act  <- act[order(act)]
  sapply(split(pred, act), tabulate, nbins=3)
}



print.plotlist<-function(xx, layout=matrix(1:length(xx)), more=F) {
  lty<-NULL
  if ( is.matrix(layout) ) {
    lyt <- layout
    col.widths <- rep.int(1, ncol(lyt))
    row.heights <- rep.int(1, nrow(lyt))
  } else if ( is.list(layout) ) {
    stopifnot(class(layout[[1]]) == "matrix")
    lyt <- layout[[1]]
    col.widths <- if (!is.null(layout$widths)) layout$widths else rep.int(1, ncol(lyt))
    row.heights <- if (!is.null(layout$heights)) layout$heights else rep.int(1, nrow(lyt))
  }
  stopifnot(length(col.widths)==ncol(lty))
  stopifnot(length(row.heights)==nrow(lty))
  maxi <- max(lyt)
  col.pts <- cumsum(c(0, col.widths))/sum(col.widths)
  row.pts <- rev(cumsum(c(0, rev(row.heights)))/sum(row.heights))
  for(i in 1:length(xx)) {
    j <-((i-1)%%maxi)+1
    wch <- which(lyt==j, arr.ind=T)
    stopifnot(nrow(wch)>0)
    pos <- apply(wch,2,range)
    ps <- c(col.pts[pos[1,2]], row.pts[pos[2,1]+1], col.pts[pos[2,2]+1],row.pts[pos[1,1]])
    print(
      xx[[i]], 
      position = ps,
      #split=c(rev(which(lyt==j, arr.ind=T)),rev(dim(lyt))),
      more=ifelse(j != maxi & i<length(xx), T, more)
    )
  }
  invisible(F)
}


GroupInt = function(listpred, X, Y, colpred="1", xlab, Ncomp)
{
  
  colpred <- match.arg(colpred)
  
  if(colpred=="1"){
    col.reg = cividis
  }
  if(colpred=="2"){
    col.reg = viridis
  }
  if(colpred=="3"){
    col.reg = magma
  }
  
  Lpred = length(listpred)
  
  
  library(lattice)
  plots<-lapply(1:Lpred, function(i) {levelplot(listpred[[i]]~X + Y,
                                                col.regions = col.reg(40),
                                                colorkey=FALSE,
                                                scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                                                xlab = xlab[i], ylab = "",              # remove axis title         # change font size for x- & y-axis text
                                                main = "")})
  
  print.plotlist(plots, layout=matrix(1:Lpred, ncol=Ncomp))
  
}

Member_prob = function(sim){
  
  nclust = length(unique(sim$Known.ppp$marks))
  
  #par(xpd=NA)
  known.marks = unique(sim$iterppp$marks)
  plot(x=seq_along(sim$Weight.df[,1]), y=sim$Weight.df[,1], col = "orange", pch=16, ylim=c(0,1),
       xlab="observations", ylab="weight")
  colvect=c("purple", "turquoise3", "darkred", "green", "brown")[1:nclust-1]
  for (i in 2:nclust) {
    points(x=seq_along(sim$Weight.df[,i]), y=sim$Weight.df[,i], col = colvect[i-1], pch=16, ylim=c(0,1))
    i =i + 1
    legend(25,0.5, c(known.marks), col = c("orange", colvect),
           pch = 16, xjust = 1, yjust = 0, merge = FALSE)
    
  }
}


Coef_fit=function(sim){
  
  if (sim$classif =="Loop"){
    classif = "soft"
  }else{
    classif = sim$classif
    
  }
  n.sp = length(unique(sim$Known.ppp$marks))
  #classif <- match.arg(classif)
  
  if(sim$classif == "hard"){
    cov.list = sim$fit.final$covariates
    Coef_mat = matrix(NA, nrow = length(cov.list)+1, ncol = n.sp)

    # Prepare / rearrange coef
    Coefall = sim$fit.final$coef
    
    Coef.rea = Coefall[-2]
    Coefout= Coefall[2]
    New.coef = c(Coef.rea[1:n.sp], Coefout, 
      Coef.rea[(n.sp+1):length(Coef.rea)])
    
    int_coef = New.coef[1:n.sp]
    other_coef = New.coef[-(1:n.sp)]
    Var_coef = other_coef[1:length(cov.list)]
    Rest_coef = other_coef[-(1:length(cov.list))]
    
    CoefRea = list()
    Coef_mat[,1] = c(int_coef[1],Var_coef)
    for (i in 1:(n.sp-1)) {
        seqnum = seq(from=i, to=length(cov.list)*(n.sp-1), by=n.sp-1)
        Coef_mat[,i+1] = c(int_coef[i+1], Rest_coef[seqnum]+Var_coef)
        
    }
    
    mat_colname = c()
    for (i in 1:n.sp) {
      mat_colname[i] = paste("Sp_", i, sep = "")
    }
    
  }
  
  if(sim$classif == "soft"){
    cov.list = sim$fit.final[[1]]$covariates
    Coef_mat = matrix(NA, nrow = length(cov.list)+1, ncol = n.sp)
    mat_colname = c()
    for (i in 1:n.sp) {
      Coef_mat[,i] = sim$fit.final[[i]]$coef
      mat_colname[i] = paste("Sp_", i, sep = "")
    }
    }
  
  coefMat = as.data.frame(Coef_mat)
  colnames(coefMat) = mat_colname
  rownames(coefMat) = c("Intercept", names(cov.list))
  coefMat
}


pred_int = function(sim, quadsenv, colpred = c("1", "2", "3")){
  
  n.sp = length(unique(sim$Known.ppp$marks))
  
  listpred = list()
  xlab = c()
  for (i in 1:n.sp) {
    listpred[[i]] = as.vector(t(simknn$pred.loc[[i]]$v))
    xlab[i] = paste("Sp_", i, sep = "")
  }
  
  
  # Choose color
  library(viridis)
  colpred <- match.arg(colpred)
  
  if(colpred=="1"){
    col.reg = cividis
  }
  if(colpred=="2"){
    col.reg = viridis
  }
  if(colpred=="3"){
    col.reg = magma
  }
  
  Lpred = length(listpred)
  library(lattice)
  plots<-lapply(1:Lpred, function(i) {levelplot(listpred[[i]]~quadsenv$X + quadsenv$Y,
                                                col.regions = col.reg(40),
                                                colorkey=FALSE,
                                                scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                                                xlab = xlab[i], ylab = "",              # remove axis title         # change font size for x- & y-axis text
                                                main = "")})
  
  print.plotlist(plots, layout=matrix(1:Lpred, ncol=n.sp, nrow=1))
  
}


se_pred = function(fit, ngrid=NULL){
  se.sim = se.range = se.simdat = list()
  seplot.sim = seplot.vec = list()  
  
  n.sp = length(unique(marks(fit$Known.ppp)))
  for (i in 1:n.sp) {
    if(fit$classif=="Loop"){
      se.sim[[i]] <- predict(fit$fit.final[[i]], locations=fit$sp_aug_ppp.list[[i]], se=TRUE)
      se.simdat[[i]] = se.sim[[i]]$se
      
      if(is.null(ngrid)){
        seplot.sim[[i]] <- predict(fit$fit.final[[i]], se=TRUE)
      }else{
        seplot.sim[[i]] <- predict(fit$fit.final[[i]], se=TRUE, ngrid=ngrid)
      }
      
      seplot.vec[[i]] = as.vector(seplot.sim[[i]]$se$v)
      se.range[[i]] = se.sim[[i]]$se
      
      Xplot = as.data.frame(fit$fitaft.pred[[1]])$x
      Yplot = as.data.frame(fit$fitaft.pred[[1]])$y
      
    }else{
      se.sim[[i]] <- predict(fit$fit.final[[i]], locations=fit$sp_aug.list[[i]], se=TRUE)
      se.simdat[[i]] = se.sim[[i]]$se
      
      if(is.null(ngrid)){
        seplot.sim[[i]] <- predict(fit$fit.final[[i]], se=TRUE)
      }else{
        seplot.sim[[i]] <- predict(fit$fit.final[[i]], se=TRUE, ngrid=ngrid)
      }

      seplot.vec[[i]] = as.vector(seplot.sim[[i]]$se$v)
      se.range[[i]] = se.sim[[i]]$se
      
      Xplot = as.data.frame(fit$fitaft.pred[[1]])$x
      Yplot = as.data.frame(fit$fitaft.pred[[1]])$y
    }
  }
  
  return(list(se.sim = se.sim, 
              se.simdat = se.simdat,
              seplot.sim = seplot.sim,
              seplot.vec = seplot.vec,
              se.range = se.range,
              Xplot = Xplot,
              Yplot= Yplot))
}

