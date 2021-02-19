library(spatstat)
library(lattice)
library(sp)
library(maptools)
library(raster)
library(geostatsp)
library(rgdal)
library(caret)
library(gridExtra)
library(grid)
library(latticeExtra)


#------------------------------------------------------------------------------
# Real data test
#------------------------------------------------------------------------------
# Myxophies data
Northcatnew <- read.csv("MyxophiesDATA.csv", header=TRUE, sep = ",")

nrow(Northcatnew)
head(Northcatnew)

Northmix.df = as.data.frame(Northcatnew)
Northmix= subset(Northmix.df, select=c(specie, X, Y))
colnames(Northmix)[2] <- "X"
colnames(Northmix)[3] <- "Y"

Northmix.new = Northmix


#-------------------------------------------------------------
# Environmental covariates
Gshape5k <- readShapeSpatial("C:/Users/Default/Desktop/OneDrive - The University Of Newcastle/ArcGIS/BandGrid5k_220118.shp")
shapeuse = Gshape5k

#5k
Gsub <- subset(shapeuse, select=c(OBJECTID_1, AUS_alt, CLIM05_1k1, CLIM06_1k, CLIM11_1k,
                                  CLIM13_1k, CLIM18_1k, NEAR_DIST, NEAR_DIS_1, X, Y))

badpts = which(Gsub$CLIM05_1k1 == 0)
badpts2 = which(Gsub$NEAR_DIST > 2000)

Gsub.df = data.frame(Gsub)  # make dataframe to use in function

Gsub.df.use1 = Gsub.df[-badpts,]
Gsub.df.use = Gsub.df.use1[-badpts2,]

Gsub.df.use$X = floor(Gsub.df.use1$X[-badpts2]/1000 + 0.001) #m 
Gsub.df.use$Y = floor(Gsub.df.use1$Y[-badpts2]/1000 + 0.001) #m

diff(sort(unique(Gsub.df.use$X)))
diff(sort(unique(Gsub.df.use$Y)))

colnames(Gsub.df.use)[1] <- "Rank"
colnames(Gsub.df.use)[8] <- "Dist_stream"
colnames(Gsub.df.use)[9] <- "Dist_road"
colnames(Gsub.df.use)[12] <- "coordX"
colnames(Gsub.df.use)[13] <- "coordY"
colnames(Gsub.df.use)[10] <- "X"
colnames(Gsub.df.use)[11] <- "Y"

# new var Gsub
Gsub.df.use$Dist_stream_sqrt = sqrt(Gsub.df.use$Dist_stream)
Gsub.df.use$Dist_road_sqrt  = sqrt(Gsub.df.use$Dist_road)
Gsub.df.use$Dist_stream_log = log(Gsub.df.use$Dist_stream)
Gsub.df.use$Dist_road_log  = log(Gsub.df.use$Dist_road)

Gsub.new = cbind(Gsub.df.use[,2:9], Gsub.df.use[,12:13], Gsub.df.use[,10:11])


##################################
##### 2. Set up study window #####
##################################

{
  ux = sort(unique(Gsub.new$X))
  uy = sort(unique(Gsub.new$Y))
  nx = length(ux)
  ny = length(uy)
  
  col.ref = match(Gsub.new$X, ux)
  row.ref = match(Gsub.new$Y, uy)
  
  Northmix.vec          = rep(NA, max(row.ref)*max(col.ref))
  vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
  Northmix.vec[vec.ref] = 1
  Northmixdf.mask      = matrix(Northmix.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
  Northmixdf.win       = as.owin(im(Northmixdf.mask, xcol = ux, yrow = uy)) # Window indicating accessible locations
  new.win = as.mask(Northmixdf.win, dimyx = c(ny, nx))
}

######################"#######################################################
# Marked ppm
##############################################################################
Gsub.c = scale(Gsub.new[,1:10], center = TRUE, scale = TRUE)
Gsubnew.c = cbind(Gsub.new[,11:12], Gsub.c)

Northmix.pppdf  = ppp(Northmix.new$X/1000, Northmix.new$Y/1000, window = new.win, marks = factor(Northmix.new$specie))

Northmix.unique.pppdf  = unique.ppp(Northmix.pppdf, marks = factor(Northmix.new$specie))
#rejects
x.reject = attributes(Northmix.pppdf)$rejects$x
y.reject = attributes(Northmix.pppdf)$rejects$y

reject.id = paste(x.reject, y.reject)
df.id     = paste(Northmix.new$X, Northmix.new$Y)
df.out    = which(df.id %in% reject.id)

Northmix.new.norejects = Northmix.new[-df.out,]


##
quads.df     = ppp(Gsubnew.c$X, Gsubnew.c$Y, window = new.win)

Qm     = quadscheme(data = Northmix.pppdf, dummy = quads.df, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))
Qmunique     = quadscheme(data = Northmix.unique.pppdf, dummy = quads.df, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))



X.des = Gsubnew.c

cov.list = list()
for (i in 1:dim(X.des)[2])
{
  Northmix.vec          = rep(NA, max(row.ref)*max(col.ref))
  vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
  Northmix.vec[vec.ref] = X.des[,i]
  v.i = im(matrix(Northmix.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux)), xcol = ux, yrow = uy)
  cov.list[[i]] = v.i
}
names(cov.list) = names(X.des)

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


##### Test before function
Mcarb.pp <- split(Northmix.pppdf)$Mcarb
Mcog.pp <- split(Northmix.pppdf)$Mcog
Msch.pp <- split(Northmix.pppdf)$Msch

test.ppp=split(Northmix.pppdf)$Unknown
marks(test.ppp)="Unknown"

true.ppp=superimpose(Mcarb=Mcarb.pp, Mcog=Mcog.pp, Msch=Msch.pp)


yrange= c(min(Mcarb.pp$y, Msch.pp$y, Mcog.pp$y, test.ppp$y), max(Mcarb.pp$y, Msch.pp$y, Mcog.pp$y, test.ppp$y))
xrange= c(min(Mcarb.pp$x, Msch.pp$x, Mcog.pp$x, test.ppp$x), max(Mcarb.pp$x, Msch.pp$x, Mcog.pp$x, test.ppp$x))

plot(Mcarb.pp$x, Mcarb.pp$y, pch=16, cex = 1.2, ylim=yrange,
     xlim=xrange, ylab="Y", xlab="X")
points(Mcarb.pp$x, Mcarb.pp$y, col="orange", pch=16, cex = 1.2)
points(Mcog.pp$x, Mcog.pp$y, col = "purple", pch=17, cex = 1.2)
points(Msch.pp$x, Msch.pp$y, col = "turquoise3", pch=15, cex = 1.2)
points(test.ppp$x, test.ppp$y, pch="?", col = "black", cex = 0.8)
legend(-220, 8270,legend=c("Mcarb", "Mcog","Msch"),
       col=c("orange", "purple", "turquoise3"), cex=1.2,
       pch=c(16,17,15), bty = "n")
legend(-220, 8190,legend=c("Unknown"),
       col=c("black"), cex=1.2,
       pch="?", bty = "n")


Gsubdf = as.data.frame(Gsub.c[,1:8])
env.xy = as.data.frame(cbind(quads.df$x, quads.df$y, Gsubdf))   #quads only needed in the context of ppp pbject implemented
colnames(env.xy) = c("X", "Y", paste0(colnames(Gsubdf)))

formS = as.formula(paste("~", paste("polynom(", names(Gsubdf)[1:7], ", 2)", collapse="+"),
                         "+", paste(names(Gsubdf)[8], collapse="+")))


#------------------------------------------------------------------------------------
#  function ppmMixEngine
#------------------------------------------------------------------------------------
#setwd("C:/Users/c3286500/Documents/SimulationProject1/Script/")
source("functionTestsim160420-SH.R")

# function to test
Testmixfunction = ppmMixEngine(Known.ppp=true.ppp, Unknown.ppp=test.ppp,
                               initweights = "knn", n.sp=3, quadsenv = env.xy,
                               k=1, ks=NULL, nstart=NULL, ppmform = formS, 
                               cov.bias=NULL, classif = "hard",
                               verbose = TRUE, tol = 0.000001, maxit = 50,
                               plots = TRUE)


Mcarb.pred = as.data.frame(Testmixfunction$fitaft.pred[[1]])
Mcog.pred = as.data.frame(Testmixfunction$fitaft.pred[[2]])
Msch.pred = as.data.frame(Testmixfunction$fitaft.pred[[3]])

#levelplot
all.val = c(Mcarb.pred$value, Mcog.pred$value, Msch.pred$value)
col.l <- colorRampPalette(c("gold", "violetred","midnightblue"))
ckey <- list(labels=list(cex=1.2))
minVal=min(Mcarb.pred$value, Mcog.pred$value, Msch.pred$value)
maxVal=max(Mcarb.pred$value, Mcog.pred$value, Msch.pred$value)
Q1Val = mean(quantile(Mcarb.pred$value)[2], quantile(Mcog.pred$value)[2],
             quantile(Msch.pred$value)[2])

Lcarb=levelplot(Mcarb.pred$value~Mcarb.pred$x + Mcarb.pred$y,
          at = unique(c(seq(minVal, 5e-2, length=10),
                       seq(5e-2, maxVal, length=30))),
          col.regions = col.l,
          colorkey=F,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          xlab = "", ylab = "",              # remove axis title         # change font size for x- & y-axis text
          main = list(label = "M. Carb",
                      cex = 1.5))

Lcog=levelplot(Mcog.pred$value~Mcog.pred$x + Mcog.pred$y,
          col.regions = col.l,
          at = unique(c(seq(minVal, 1e-2, length=10),
                         seq(1e-2, maxVal, length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=F,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. Cog",
                      cex = 1.5))

Lsch=levelplot(Msch.pred$value~Msch.pred$x + Msch.pred$y,
          col.regions = col.l,
          at = unique(c(seq(minVal, 1e-2, length=10),
                        seq(1e-2,maxVal, length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=ckey,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. sch",
                      cex = 1.5))


comb_levObj <- c(Lcarb, Lcog, Lsch, 
                 layout = c(3, 1), merge.legends = T)

update(comb_levObj,
       xlab = c("Mcarb", "Mcog", "Msch"),
       main="Mixture models - knn")

minval1=minVal
maxval1=maxVal

### Individual model
#-----------------------------------------------------------------
# M carb
Mcarb.pp <- split(Northmix.pppdf)$Mcarb
unmark(Mcarb.pp)

##
Qm.Mcarb = quadscheme(data = Mcarb.pp, dummy = quads.df, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))
ft.Mcarb  = ppm(Qm.Mcarb, trend = formS, covariates = cov.list, gcontrol = list(epsilon = 1e-6, maxit = 100))

DMcarb = diagnose.ppm(ft.Mcarb, type = "Pearson", main="diagnose - ft.Mcarb")
Mcarb.pred2b  = predict(ft.Mcarb, covariates = cov.list)
par(mfrow=c(1,1))
plot(log(Mcarb.pred2b), main="predict - log(ft.Mcarb2b)")
plot(Mcarb.pred2b, main="predict - ft.Mcarb2b")

#-----------------------------------------------------------------
# M cog
Mcog.pp <- split(Northmix.pppdf)$Mcog
unmark(Mcog.pp)

##
Qm.Mcog = quadscheme(data = Mcog.pp, dummy = quads.df, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))
ft.Mcog  = ppm(Qm.Mcog, trend = formS, covariates = cov.list, gcontrol = list(epsilon = 1e-6, maxit = 100))

DMcog = diagnose.ppm(ft.Mcog, type = "Pearson", main="diagnose - ft.Mcog2b")
Mcog.pred2b  = predict(ft.Mcog, covariates = cov.list)
plot(log(Mcog.pred2b), main="predict - log(ft.Mcog2b)")
plot(Mcog.pred2b, main="predict - ft.Mcog2b")

#--------------------------------------------------------------------
# M sch
Msch.pp <- split(Northmix.pppdf)$Msch
unmark(Msch.pp)

##
Qm.Msch = quadscheme(data = Msch.pp, dummy = quads.df, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))
ft.Msch  = ppm(Qm.Msch, trend = formS, covariates = cov.list, gcontrol = list(epsilon = 1e-6, maxit = 100))

DMsch = diagnose.ppm(ft.Msch, type = "Pearson", main="diagnose - ft.Msch2b")
Msch.pred2b  = predict(ft.Msch, covariates = cov.list)
plot(log(Msch.pred2b), main="predict - log(ft.Msch2b)")
plot(Msch.pred2b, main="predict - ft.Msch2b")

Mcarbdf = as.data.frame(Mcarb.pred2b)
Mcogdf = as.data.frame(Mcog.pred2b)
Mschdf = as.data.frame(Msch.pred2b)

all.valind = c(Mcarbdf$value, Mcogdf$value, Mschdf$value)
col.l <- colorRampPalette(c("gold", "violetred","midnightblue"))
ckey <- list(labels=list(cex=1.2))
minValind=min(Mcarbdf$value, Mcogdf$value, Mschdf$value)
maxValind=max(Mcarbdf$value, Mcogdf$value, Mschdf$value)


L1ind=levelplot(Mcarbdf$value~Mcarbdf$x + Mcarbdf$y,
          col.regions = col.l,
          at = unique(c(seq(min(all.valind), 5e-3, length=10),
                        seq(5e-3, max(all.valind), length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=F,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. Carb",
                      cex = 1.5))

L2ind=levelplot(Mcogdf$value~Mcogdf$x + Mcogdf$y,
          col.regions = col.l,
          at = unique(c(seq(min(all.valind), 5e-3, length=10),
                        seq(5e-3, max(all.valind), length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=F,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. Cog",
                      cex = 1.5))

L3ind=levelplot(Mschdf$value~Mschdf$x + Mschdf$y,
          col.regions = col.l,
          at = unique(c(seq(min(all.valind), 5e-3, length=10),
                        seq(5e-3, max(all.valind), length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=ckey,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. sch",
                      cex = 1.5))


comb_levObj <- c(L1ind, L2ind, L3ind, 
                 layout = c(3, 1), merge.legends = T)

update(comb_levObj,
       xlab = c("Mcarb", "Mcog", "Msch"),
       main="Individual models")

cor12=cor(Mcarbdf$value, Mcogdf$value)
cor13=cor(Mcarbdf$value, Mschdf$value)
cor23=cor(Mcogdf$value, Mschdf$value)


#-----------------------------------------------------------------------------------------------------------------------
quads.win = new.win

Testrealadd = ppmLoopEngine(Known.ppp=true.ppp, Unknown.ppp=test.ppp,
                            n.sp=3, quadsenv = env.xy,
                            addpt = "LoopT", ppmform=formS,
                            delta_max=0.5, delta_min=0.1, delta_step =0.1,
                            num.add = NULL, cov.bias=NULL, SetValBias =NULL, 
                            rAreaInter=NULL, maxit = 50,
                            tol=0.000001, verbose = TRUE, plots = FALSE)

par(xpd=T, mar=par()$mar+c(0,0,1,0))
pl.weights = Testrealadd$New_weights
plot(x=seq_along(pl.weights[,1]), y=pl.weights[,1], col = "orange", pch=16, ylim=c(0,1),
     xlab="observations", ylab="weight")
points(x=seq_along(pl.weights[,2]), y=pl.weights[,2], col = "purple", pch=18, ylim=c(0,1))
points(x=seq_along(pl.weights[,3]), y=pl.weights[,3], col = "Turquoise3", pch=17, ylim=c(0,1))
legend(1,1, c("Mcarb", "Mcog", "Msch"), col = c("orange", "purple", "Turquoise3"),
       pch = c(16, 18, 17), xjust = 1, yjust = 0, merge = FALSE)



# Loop method

McarbL.df  = as.data.frame(predict(Testrealadd$fit.final[[1]]))
McogL.df  = as.data.frame(predict(Testrealadd$fit.final[[2]]))
MschL.df  = as.data.frame(predict(Testrealadd$fit.final[[3]]))


all.valL = rbind(McarbL.df$value, McogL.df$value, MschL.df$value)
minVal=min(McarbL.df$value, McogL.df$value, MschL.df$value)
maxVal=max(McarbL.df$value, McogL.df$value, MschL.df$value)
#Q1Val=quantile(all.val)[2]
col.l <- colorRampPalette(c("gold", "violetred","midnightblue"))
ckey <- list(labels=list(cex=1.2))


#other scale
L1loop=levelplot(McarbL.df$value~McarbL.df$x + McarbL.df$y,
          col.regions = col.l,
          at = unique(c(seq(minVal, 1e-3, length=10),
                        seq(1e-3, maxVal, length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=F,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. Carb",
                      cex = 1.5))

L2loop=levelplot(McogL.df$value~McogL.df$x + McogL.df$y,
          col.regions = col.l,
          at = unique(c(seq(minVal, 1e-3, length=10),
                        seq(1e-3, maxVal, length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=F,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. Cog",
                      cex = 1.5))

L3loop=levelplot(MschL.df$value~MschL.df$x + MschL.df$y,
          col.regions = col.l,
          at = unique(c(seq(minVal, 1e-3, length=10),
                        seq(1e-3, maxVal, length=30))),
          xlab = "", ylab = "",              # remove axis titles
          colorkey=ckey,
          scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
          main = list(label = "M. sch",
                      cex = 1.5))


comb_levObj <- c(L1loop, L2loop, L3loop, 
                 layout = c(3, 1), merge.legends = T)
#print(comb_levObj)
update(comb_levObj,
       xlab = c("Mcarb", "Mcog", "Msch"),
       main="Loopgr models")


## all plots

minVal=min(minval1, minValind, minVal)
maxVal=max(maxval1, maxValind, maxVal)

Lcarb=levelplot(Mcarb.pred$value~Mcarb.pred$x + Mcarb.pred$y,
                at = unique(c(seq(minVal, 1e-3, length=10),
                              seq(1e-3, maxVal, length=15))),
                col.regions = col.l,
                colorkey=F,
                scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                xlab = "", ylab = "",              # remove axis title         # change font size for x- & y-axis text
                main = list(label = "M. Carb",
                            cex = 1.5))

Lcog=levelplot(Mcog.pred$value~Mcog.pred$x + Mcog.pred$y,
               col.regions = col.l,
               at = unique(c(seq(minVal, 1e-3, length=10),
                             seq(1e-3, maxVal, length=15))),
               xlab = "", ylab = "",              # remove axis titles
               colorkey=F,
               scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
               main = list(label = "M. Cog",
                           cex = 1.5))

Lsch=levelplot(Msch.pred$value~Msch.pred$x + Msch.pred$y,
               col.regions = col.l,
               at = unique(c(seq(minVal, 1e-3, length=10),
                             seq(1e-3,maxVal, length=15))),
               xlab = "", ylab = "",              # remove axis titles
               colorkey=F,
               scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
               main = list(label = "M. sch",
                           cex = 1.5))


L1ind=levelplot(Mcarbdf$value~Mcarbdf$x + Mcarbdf$y,
                col.regions = col.l,
                at = unique(c(seq(minVal, 1e-3, length=10),
                              seq(1e-3, maxVal, length=15))), 
                xlab = "", ylab = "",# remove axis titles
                colorkey=F,
                scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                main = list(label = "M. Carb",
                            cex = 1.5))

L2ind= levelplot(Mcogdf$value~Mcogdf$x + Mcogdf$y,
                 col.regions = col.l,
                 at = unique(c(seq(minVal, 1e-3, length=10),
                               seq(1e-3, maxVal, length=15))), 
                 xlab = "", ylab = "",# remove axis titles
                 colorkey=F,
                 scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                 main = list(label = "M. Cog",
                             cex = 1.5))

L3ind=levelplot(Mschdf$value~Mschdf$x + Mschdf$y,
                col.regions = col.l,
                at = unique(c(seq(minVal, 1e-3, length=10),
                              seq(1e-3, maxVal, length=15))), 
                xlab = "", ylab = "",              # remove axis titles
                colorkey=F,
                scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                main = list(label = "M. sch",
                            cex = 1.5))

L1loop=levelplot(McarbL.df$value~McarbL.df$x + McarbL.df$y,
                 col.regions = col.l,
                 at = unique(c(seq(minVal, 1e-3, length=10),
                               seq(1e-3, maxVal, length=15))),
                 xlab = "", ylab = "",              # remove axis titles
                 colorkey=F,
                 scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                 main = list(label = "M. Carb",
                             cex = 1.5))

L2loop=levelplot(McogL.df$value~McogL.df$x + McogL.df$y,
                 col.regions = col.l,
                 at = unique(c(seq(minVal, 1e-3, length=10),
                               seq(1e-3, maxVal, length=15))),
                 xlab = "", ylab = "",              # remove axis titles
                 colorkey=F,
                 scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                 main = list(label = "M. Cog",
                             cex = 1.5))

L3loop=levelplot(MschL.df$value~MschL.df$x + MschL.df$y,
                 col.regions = col.l,
                 at = unique(c(seq(minVal, 1e-3, length=10),
                               seq(1e-3, maxVal, length=15))),
                 xlab = "", ylab = "",              # remove axis titles
                 colorkey=ckey,
                 scales = list(y = list(draw = FALSE), x = list(draw = FALSE)),
                 main = list(label = "M. sch",
                             cex = 1.5))


comb_levObj <- c(L1loop, L2loop, L3loop, 
                 L1ind, L2ind, L3ind,
                 Lcarb, Lcog, Lsch, 
                 layout = c(3, 3), merge.legends = T)
#print(comb_levObj)
update(comb_levObj,
       xlab = c("Mcarb", "Mcog", "Msch"),
       ylab = c("LoopT", "indiv", "Mixture knn"),
       main="Myxophies species predicted distribution")
