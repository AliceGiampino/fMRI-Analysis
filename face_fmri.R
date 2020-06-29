
rm(list = ls())

# setwd("path")

# import data face_rfx > cons_fir
# Subject-specific t-contrasts 
# on the main effect of faces versus baseline
# on the canonical HRF

source("threshold_k.R")
source("optimal_cluster.R")
source("cluster_gini.R")
source("threshold.R")


library(fmri)
PATH = 'cons_can/'
ds=read.ANALYZE(paste(PATH,"con_0", sep=""), 
                numbered=TRUE,
                postfix="",
                picstart=6,
                numbpic=12)


# alternative: face_rfx > cons_fir
# for the main effect vs baseline 

# another alternative: face_rfx > cons_informed
# The contrast images for the difference 
# between famous and nonfamous faces happen 
# to correspond to contrast numbers 39-50 
# in the 1st-level model


library(AnalyzeFMRI)
allfiles = dir('cons_can', 
               pattern = '.hdr', 
               full.names = TRUE)
files = allfiles
imageDim = f.read.analyze.header(files[1])$dim[2:4]
mask3D = array(1, imageDim)
for (file in files) {
  img =f.read.analyze.volume(file)[,,,1]
  mask3D = mask3D*(!is.na(img))
}
mask_ind = which(mask3D == 1)

#============
# BRODMANN
#============

mask.brodmann=f.read.nifti.volume("E:/PhD/Statistical Inference II/HW5/brodmann_my.nii")
mask2 = (mask3D==1 & mask.brodmann[,,,1]>0)

mask_ind2 = which(mask2 == 1)

table(mask3D, mask.brodmann>0)

require("rgl")
require("misc3d")
# zoom = 1
# windowRect = c(0,45,630,630)
# open3d(zoom = zoom, windowRect=windowRect)
# rgl.viewpoint( theta = 0, phi = 0)
# contour3d(mask3D, level=1, alpha=0.1)
# map.broad <- (mask.brodmann[,,,1]>0)==T
# map.broad[mask3D==0] <- 0
# map.broad[map.broad==F] <- 0
# map.broad[map.broad==T] <- 1
# contour3d(map.broad, level=1, alpha=0.7, add=T, color="red")
# dmap <- dim(mask3D)
# text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
# text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
# text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
# text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")

#=============
# data matrix 
#=============

ds4D = extractData(ds)
X = NULL
n = ds$dim[4]
for(i in 1:n) {
  scan = ds4D[,,,i]
  X = rbind(X, scan[mask_ind])
}

# one of the scan:
scan <- ds4D[,,,1]
x1 <- scan[mask_ind]
zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, level=1, alpha=0.1)
contour3d(scan, level=1, add=T, color="red", alpha=0.8)
dmap <- dim(mask3D)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


#=============
# ARI
#=============

library(hommel)
thetahat3D = apply(ds4D, c(1,2,3), function(x)
  mean(x)/sd(x)
)

thetahat = apply(X, 2, function(x) {
  mean(x)/sd(x) })
m = length(thetahat)

tmap = apply(X, 2, function(x) {
  t.test(x)$statistic })

pval = apply(X, 2, function(x) {
  t.test(x)$p.value }) 

res=hommel(pval) # adjusted pvalues of hommel

summary(res, alpha=0.05)
# With 0.95 confidence: at least 7745 discoveries.
# 248 hypotheses with adjusted p-values below 0.05.

# retrieve familywise error adjusted p-values
p.adj <- hommel::p.adjust(res)

# Find the concentration set bound
threshold <- concentration(res) 
# The concentration set is the subset of the p-values that
# contains all discoveries at confidence level 1-alpha
# Returns a p-value. P-values larger than that value contain no discoveries 
# in this data set at this level of alpha and may be disregarded.

# Find the concentration set itself
set <- pval <= threshold
sum(set)


localtest(res, tdp=threshold)
# This function calculates familywise error adjusted
# p-values for local tests. Familywise error control 
# is over all possible local tests simultaneously

# P-value
sum(pval <= 0.05) # 21870
bonf <- p.adjust(pval, method="bonferroni")
pval_B <- bonf <= 0.05
sum(pval_B) # 212
pval_hommel <-  res@adjusted <= 0.05
sum(pval_hommel) # 248
BH <- p.adjust(pval, method="BH")
pval_BH <- BH <= 0.05
sum(pval_BH) # 16258
BY <- p.adjust(pval, method="BY")
pval_BY <- BY <= 0.05
sum(pval_BY) # 5813

# Visualization of one scan:
scan <- ds4D[,,,1]
scan <- thetahat3D
scan[scan=="NaN"] <- 0
x1 <- scan[mask_ind]
contour3d(x1, level=1, add=T, color="blue")


# visualize pvalues <= 0.05 into the brain
sum(scan==0)
scan_hom <- scan
scan_hom[scan!=0] <- pval_hommel 
scan_BH <- scan
scan_B <- scan
scan_BY <- scan
scan_BH[scan!=0] <- pval_BH
scan_B[scan!=0] <- pval_B
scan_BY[scan!=0] <- pval_BY


# Bonferroni and Hommel
zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, level=1, alpha=0.1)
contour3d(scan_hom, level=1, add=T, color="red", alpha=0.8) # hommel
contour3d(scan_B, level=1, add=T, color="blue") # bonferroni
dmap <- dim(mask3D)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


# BH and BY
zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, level=1, alpha=0.1)
contour3d(scan_BH, level=1, add=T, color="brown1", alpha=0.5)
contour3d(scan_BY, level=1, add=T, color="yellow", alpha=1)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


#=============
# CLUSTER
#=============


nmat <- expand.grid(-1:1, -1:1, -1:1)
nmat <- nmat[-c(1,3,7,9,14,19,21,25,27), ]
stat = n*(thetahat3D)^2
stat[mask3D==F] = -1

clusters = cluster.threshold(stat, nmat=nmat, 
                             level.thr = 30, # first try for level threshold
                             size.thr = 50)
clusters= clusters[mask3D==1]==1

# higher the threshold, lower the cluster number

coord=which(mask3D !=0, arr.ind = TRUE) # coordinates where the mask is not 0 (the brain)
x = coord[,1]
y = coord[,2]
z = coord[,3]

dd = dist(coord[clusters,])
hc = hclust(dd, "single")
plot(hc)
abline(h=1, col=2)
K=6
ct = cutree(hc,k=K)
table(ct)

clusters2 = clusters
clusters2[clusters] = ct


map <- mask3D
for (i in 1:nrow(coord)){
  map[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
}


for (j in 1:K){
  id_selected=which( map[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", discoveries(res,alpha=.05, ix=id_selected), "(", discoveries(res,alpha=0.01, ix=id_selected), ")"), "\n")
  cat(paste("%:", round(tdp(res,alpha=.05, ix=id_selected)*100,1), "%", "(", round(tdp(res,alpha=0.01, ix=id_selected)*100,1), "% )"), "\n")
}

# number of discoveries (tdp = true discoveries proportion), a lower bound
# ( fdp upper bound )

# 2561 
# N: 2393 ( 1682 ) 
# %: 93.4 % ( 65.7 % ) 
# 643 
# N: 477 ( 102 ) 
# %: 74.2 % ( 15.9 % ) 
# 83 
# N: 0 ( 0 ) 
# %: 0 % ( 0 % ) 
# 78 
# N: 16 ( 5 ) 
# %: 20.5 % ( 6.4 % ) 
# 347 
# N: 195 ( 47 ) 
# %: 56.2 % ( 13.5 % ) 
# 172 
# N: 23 ( 0 ) 
# %: 13.4 % ( 0 % ) 

# Bonferroni attivation ----

for (j in 1:K){
  id_selected=which( map[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", sum(pval_B[id_selected])), "\n")
  cat(paste("%:", round(sum(pval_B[id_selected])/length(id_selected)*100,2)), "\n")
}
# 2561 
# N: 183 
# %: 7.15 
# 643 
# N: 7 
# %: 1.09 
# 83 
# N: 0 
# %: 0 
# 78 
# N: 8 
# %: 10.26 
# 347 
# N: 14 
# %: 4.03 
# 172 
# N: 0 
# %: 0 

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
colsK =  rainbow(K)
for (i in 1:length(colsK)){
  contour3d( map==i, alpha = 0.5, add = TRUE, draw=TRUE, level=0.5, color=colsK[i])
}
dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")

# Clusters ----

clusters = cluster.threshold(stat, nmat=nmat, 
                             level.thr = 45.53144, 
                             size.thr = 50)
clusters= clusters[mask3D==1]==1

coord=which(mask3D !=0, arr.ind = TRUE) 
x = coord[,1]
y = coord[,2]
z = coord[,3]

dd = dist(coord[clusters,])
hc = hclust(dd, "single")
plot(hc)
abline(h=3, col=2)
K=4
ct = cutree(hc,k=K)
table(ct)

clusters2 = clusters
clusters2[clusters] = ct


map <- mask3D
for (i in 1:nrow(coord)){
  map[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
}


for (j in 1:K){
  id_selected=which( map[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", discoveries(res,alpha=.05, ix=id_selected), "(", discoveries(res,alpha=0.01, ix=id_selected), ")"), "\n")
  cat(paste("%:", round(tdp(res,alpha=.05, ix=id_selected)*100,1), "%", "(", round(tdp(res,alpha=0.01, ix=id_selected)*100,1), "% )"), "\n")
}

# number of discoveries (tdp = true discoveries proportion) sono lower bound
# fdp upper bound, non riportato

# cluster number of voxel, alpha=0.05 and (alpha = 0.01)
# 1409 
# N: 1382 ( 1258 ) 
# %: 98.1 % ( 89.3 % ) 
# 93 
# N: 66 ( 7 ) 
# %: 71 % ( 7.5 % ) 
# 89 
# N: 63 ( 11 ) 
# %: 70.8 % ( 12.4 % ) 
# 172 
# N: 146 ( 47 ) 
# %: 84.9 % ( 27.3 % ) 

# Bonferroni attivation:

for (j in 1:K){
  id_selected=which( map[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", sum(pval_B[id_selected])), "\n")
  cat(paste("%:", round(sum(pval_B[id_selected])/length(id_selected)*100,2)), "\n")
  print(summary(pval[id_selected]))
}
# 1409 
# N: 183 
# %: 12.99 
# 93 
# N: 1 
# %: 1.08 
# 89 
# N: 5 
# %: 5.62 
# 172 
# N: 14 
# %: 8.14 

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
colsK =  rainbow(K)
for (i in 1:length(colsK)){
  contour3d( map==i, alpha = 0.5, add = TRUE, draw=TRUE, level=0.5, color=colsK[i])
}

dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")

# Brodmann areas in the selected clusters ----

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, level=1, alpha=0.1)
map.broad_sel <- (mask.brodmann[,,,1]>0)==T
map.broad_sel[mask3D==0] <- 0
map.broad_sel[map.broad==F] <- 0
map.broad_sel[map.broad==T] <- 1
map.broad_sel[map_sub_att==0] <- 0
contour3d(map.broad_sel, level=1, alpha=1, add=T, color="blue")
for (i in 1:length(colsK)){
  contour3d( map_sub_att==i, alpha = 0.3, add = TRUE, draw=TRUE, level=0.5, color=colsK[i])
}
dmap <- dim(mask3D)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


# Sub-cluster Gini ----

library(DescTools)
Gini(stat[-which(stat==-1)]) # 0.6857752

# opt <- optimal_cluster(map, coord, hc, 0.1, 4, 0.01, clusters, res)

library(dplyr)
# saveRDS(opt, "opt.Rds")
opt <- readRDS("opt.Rds")

gini_ind <- as.data.frame(sapply(opt, function(i) i$gini))
colnames(gini_ind) <- "gini"
table(gini_ind$gini)
# first 90 entries: max Gini

h <- sapply(opt, function(i) i$h) 
attivation <- sapply(opt, function(i) i$attivation)

att <- t(as.data.frame(attivation[91]))
colnames(att) <- c("discoveries", "tdp")
table(att[,1]) # discoveries
table(att[,2]) # tdp


# Choose h = 1
h[91]
gini_ind[91,]
ct_sub = cutree(hc, h=h[91])

table(ct_sub)

clusters_sub = clusters
clusters_sub[clusters] = ct_sub


map_sub <- mask3D
for (i in 1:nrow(coord)){
  map_sub[coord[i,1],coord[i,2],coord[i,3]] <- clusters_sub[i]
}

K_sub <- length(unique(ct_sub))

for (j in 1:K_sub){
  id_selected=which( map_sub[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", discoveries(res,alpha=.05, ix=id_selected)), "\n")
  cat(paste("%:", round(tdp(res,alpha=.05, ix=id_selected)*100,1), "%"), "\n")
}

# 1383 
# N: 1356 
# %: 98 % 
# 2 
# N: 0 
# %: 0 % 
#   1 
# N: 0 
# %: 0 % 
#   12 
# N: 0 
# %: 0 % 
#   92 
# N: 65 
# %: 70.7 % 
#   6 
# N: 0 
# %: 0 % 
#   4 
# N: 0 
# %: 0 % 
#   89 
# N: 63 
# %: 70.8 % 
#   1 
# N: 0 
# %: 0 % 
#   1 
# N: 0 
# %: 0 % 
#   172 
# N: 146 
# %: 84.9 %

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
colsK_sub =  rainbow(K_sub)
for (i in 1:length(colsK_sub)){
  contour3d( map_sub==i, alpha = 0.5, add = TRUE, draw=TRUE, level=0.5, color=colsK_sub[i])
}
dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
map_sub_att <- map_sub
map_sub_att[map_sub==2] <- 0
map_sub_att[map_sub==3] <- 0
map_sub_att[map_sub==4] <- 0
map_sub_att[map_sub==6] <- 0
map_sub_att[map_sub==7] <- 0
map_sub_att[map_sub==9] <- 0
map_sub_att[map_sub==10] <- 0
table(map_sub_att)

map_sub_att[map_sub==5] <- 2
map_sub_att[map_sub==8] <- 3
map_sub_att[map_sub==11] <- 4

for (i in 1:length(colsK)){
  contour3d( map_sub_att==i, alpha = 0.5, add = TRUE, draw=TRUE, level=0.5, color=colsK[i])
}
dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


# optimization of the cluster threshold -------------------------------------------------

# th <- threshold(stat, 40, 80, 1, 50, coord, 4)

thresh <- as.numeric(lapply(th, function(i) i$mean_tdp))
max(thresh)

thr <- seq(40, 80, 1)
thr[which.max(thresh)]

clusters = cluster.threshold(stat, nmat=nmat, 
                             level.thr = 69, 
                             size.thr = 50)
clusters= clusters[mask3D==1]==1

coord=which(mask3D !=0, arr.ind = TRUE) 
x = coord[,1]
y = coord[,2]
z = coord[,3]

dd = dist(coord[clusters,])
hc = hclust(dd, "single")
# plot(hc)
# abline(h=3, col=2)
K=4
ct = cutree(hc,k=K)
table(ct)

clusters2 = clusters
clusters2[clusters] = ct


map <- mask3D
for (i in 1:nrow(coord)){
  map[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
}


for (j in 1:K){
  id_selected=which( map[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", discoveries(res,alpha=.05, ix=id_selected), "(", discoveries(res,alpha=0.01, ix=id_selected), ")"), "\n")
  cat(paste("%:", round(tdp(res,alpha=.05, ix=id_selected)*100,1), "%", "(", round(tdp(res,alpha=0.01, ix=id_selected)*100,1), "% )"), "\n")
}
# 282 
# N: 279 ( 261 ) 
# %: 98.9 % ( 92.6 % ) 
# 142 
# N: 139 ( 121 ) 
# %: 97.9 % ( 85.2 % ) 
# 90 
# N: 87 ( 69 ) 
# %: 96.7 % ( 76.7 % ) 
# 57 
# N: 54 ( 36 ) 
# %: 94.7 % ( 63.2 % ) 

perc <- c()
for (j in 1:K){
  id_selected=which( map[mask3D==1]==j)
  perc[j] <- round(tdp(res,alpha=.05, ix=id_selected)*100,1)
  print(summary(pval[id_selected]))
}

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
colsK =  rainbow(K)
for (i in 1:length(colsK)){
  contour3d( map==i, alpha = 0.5, add = TRUE, draw=TRUE, level=0.5, color=colsK[i])
}
dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


# Adjusted pvalues in the clusters ----

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, level=1, alpha=0.1)
map_pval <- mask3D
for (i in 1:length(pval_hommel)){
  if(pval_hommel[i]==T){
    map_pval[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
  }
  else{
    map_pval[coord[i,1],coord[i,2],coord[i,3]] <- 0
  }
}
for (i in 1:length(colsK)){
  contour3d( map==i, alpha = 0.3, add = TRUE, draw=TRUE, level=0.5, color=colsK[i])
}
for (i in 1:length(colsK)){
  contour3d( map_pval==i, alpha = 1, add = TRUE, draw=TRUE, level=0.5, color="blue")
}

text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")

# activation by hommel:

for(j in 1:K){
  
  id_selected=which( map[mask3D==1]==j)
  
  cat(length(id_selected),"\n")
  
  pval_cluster <- res@adjusted[id_selected]
  
  m_cluster <- length(pval_cluster)
  
  cat(paste(round(conc_pvalue(pval_cluster, m_cluster)*100,1), "%"), "\n")
  
}

for(j in 1:K){
  
  id_selected=which( map[mask3D==1]==j)
  
  cat(paste("N: ", length(id_selected),"\n"))
  
  pval_cluster <- res@adjusted[id_selected] <= 0.05
  
  m_cluster <- length(pval_cluster)
  
  cat(paste("conc: ",round(conc_pvalue(pval_cluster, m_cluster)*100,1), "%"), "\n")
  
  cat(paste("tdp: ",round(tdp(res, alpha=0.05, ix=id_selected)*100,1), "%"), "\n")
  
  
}

# Bonferroni:

for(j in 1:K){
  
  id_selected=which( map[mask3D==1]==j)
  
  cat(length(id_selected),"\n")
  
  pval_cluster <- pval_B[id_selected]
  
  m_cluster <- length(pval_cluster)
  
  cat(paste(round(conc_pvalue(pval_cluster, m_cluster)*100,1), "%"), "\n")
  
}

# Calculate Gini in the clusters:
# gini <- Cluster_Gini(hc, min_h = 0.1, max_h=6, by=0.1, clusters)

gini_out <- as.numeric(lapply(gini, function(i) i$gini_out))
table(gini_out)
which.max(gini_out) 

h <- seq(0.1, 6, by=0.1)
gini[[10]]
h[10]

# considering only pvalues ----

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
colsK =  rainbow(K)
p.hom50 <- p.adj < 0.5
p.hom30 <- p.adj < 0.3
p.hom10 <- p.adj < 0.1
map_pval_hom50 <- mask3D
map_pval_hom50[mask3D!=0] <- p.hom50
map_pval_hom30 <- mask3D
map_pval_hom30[mask3D!=0] <- p.hom30
map_pval_hom10 <- mask3D
map_pval_hom10[mask3D!=0] <- p.hom10

cluster.p <- rep(0, m)
cluster.p[pval_hommel] <- 1
mapp <- mask3D
for (i in 1:nrow(coord)){
  mapp[coord[i,1],coord[i,2],coord[i,3]] <- cluster.p[i]
}

id_selected=which( mapp[mask3D==1]==1)
percent <- tdp(res, ix=id_selected, alpha=0.05)

contour3d(map_pval_hom50, alpha = 0.3, draw = TRUE, level=1, add=T, color="grey")
contour3d(scan_hom, alpha = 1.2, draw = TRUE, level=1, add=T, color=2)
dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")


# optimization of h and k ----------------------------------------------------

library(tictoc)

tic()
thk <- threshold_k(stat, 40, 90, 0.5, 50, coord, c(4:7))
toc()

thresh_h <- lapply(thk, function(i) lapply(i, function(j) lapply(j, function(w) w$mean_tdp)))
thresh_k <- lapply(thk, function(i) lapply(i, function(j) lapply(j, function(w) w$k)))

thr <- seq(40, 90, 0.5)

x <- sapply(thresh_h, function(i) i$mean_tdp)
df_h <- data.frame(matrix(unlist(x), nrow=4))

h <- thr[37]
clusters_t = cluster.threshold(stat, nmat=nmat, 
                               level.thr = h, 
                               size.thr = 100)
clusters_t = clusters_t[mask3D==1]==1

dd_t = dist(coord[clusters_t,])
hc = hclust(dd_t, "single")

plot(hc)
abline(h=1.5, col=2)
ct = cutree(hc, h=1.3)
table(ct)

clusters2 = clusters_t
clusters2[clusters_t] = ct


map <- mask3D
for (i in 1:nrow(coord)){
  map[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
}


for (j in 1:length(unique(ct))){
  id_selected=which( map[mask3D==1]==j)
  cat(length(id_selected),"\n")
  cat(paste("N:", discoveries(res,alpha=.05, ix=id_selected), "(", discoveries(res,alpha=0.01, ix=id_selected), ")"), "\n")
  cat(paste("%:", round(tdp(res,alpha=.05, ix=id_selected)*100,1), "%", "(", round(tdp(res,alpha=0.01, ix=id_selected)*100,1), "% )"), "\n")
}
# 467 
# N: 459 ( 418 ) 
# %: 98.3 % ( 89.5 % ) 
# 225 
# N: 217 ( 176 ) 
# %: 96.4 % ( 78.2 % ) 
# 163 
# N: 155 ( 115 ) 
# %: 95.1 % ( 70.6 % ) 

zoom = 1
windowRect = c(0,45,630,630)
open3d(zoom = zoom, windowRect=windowRect)
rgl.viewpoint( theta = 0, phi = 0)
contour3d(mask3D, alpha = 0.1, draw = TRUE, level=1)
colsK =  rainbow(K)
for (i in 1:length(colsK)){
  contour3d( map==i, alpha = 0.5, add = TRUE, draw=TRUE, level=0.5, color = colsK[i])
}
dmap = dim(map)
text3d(x=dmap[1]*0.99, y=dmap[2]/2, z = dmap[3]/2, text="R")
text3d(x=dmap[1]*0.01, y=dmap[2]/2, z = dmap[3]/2, text="L")
text3d(x=dmap[1]/2, y=dmap[2]*0.99, z = dmap[3]/2, text="A")
text3d(x=dmap[1]/2, y=dmap[2]*0.01, z = dmap[3]/2, text="P")

