
rm(list=ls())

setwd("E:\\PhD\\Statistical Inference II\\HW5\\face_rfx")

library(fmri)
library(AnalyzeFMRI)
library(RNifti)
library(hommel)
library(plyr)

brodmann <- read.NIFTI("E:/PhD/Statistical Inference II/HW5/brodmann_my.nii")

file <- read.ANALYZE("cons_can/con_0", numbered = T, postfix="", 6, 12)

# create expected BOLD signal and design matrix
hrf <- fmri.stimulus(scans=12, onsets=mean(0.2, 4), durations=c(0.5)); hrf

dhrf <- (c(0,diff(hrf)) + c(diff(hrf),0))/2

x <- fmri.design(cbind(hrf, dhrf))

mask <- file$mask

spm <- fmri.lm(file, x, mask=mask) # problems

# structure adaptive smoothing with maximum bandwith hmax
spmsmooth <- fmri.smooth(spm)

# calculate p-values for smoothed parametric map
pvalue <- fmri.pvalue(spmsmooth)

# plot result slicewise into file or ...
plot(pvalue, maxpvalue=0.05, device="jpeg", file="result.jpeg")

# ... plot interactive 3D visualization
plot(pvalue, maxpvalue=0.05, type="3d")


# Others ----

ttt <- extractData(test)

mean(ttt)
sd(ttt)
hist(ttt)
summary(ttt)
image(as.nifti(ttt))

# prova 1 (NON VA IL MODELLO):
img <- file
onsets <- c(5.25, 21.45,
            100.12, 1)
dur <-  c(5, 5, 10, 1)
hrf <-  fmri.stimulus(scans = file$dim0[4],
                      times = onsets,
                      durations = c(1))

x <-  fmri.design(hrf, order =2)
mask <-  file$mask
dataMatrix <-  NULL
noScans <-  test$dim[4]
for(t in 1:noScans){
  scan <-  file[,,,t]
  dataMatrix <-  rbind(dataMatrix, scan[mask])
}
lms <-  apply(dataMatrix, 2, function(y) lm(y,x))
# obtain the summary statistics for each voxel
summaries <-  lapply(lms, summary)
# Coefficients
coefs <-  lapply(summaries, function(x) x$coefficients)
# Intercept map
int_vals <-  sapply(coefs, function(x) x["(Intercept)", ])
int_est <-  int_vals["Estimate",]
int_pvals <-  int_vals["Pr(.|t|)",]
# Slope map
beta_vals <-  sapply( coefs, function(x) x[2,])
beta_est <-  beta_vals["Estimate",]
beta_pvals <-  beta_vals["Pr(.|t|)",]
# apply the function t.test to every column (voxel)
ttest.out <-  apply(dataMatrix, 2,
                    function(x) {
                      temp <-  t.test(x)
                      c(temp$statistic, temp$p.value)})
# label the rows
rownames(ttest.out) <-  c("statistic", "pvalue")

Threshold <-  Threshold.FDR(ttest.out[1,],
                            q =.05,
                            cV.type = 2,
                            type = "t",
                            df1 = noScans - 1)

outputImage <-  array(0, imageDim)
toSave <-  ttest.out[1,] >= Threshold
outputImage[mask[toSave]] <- ttest.out[1, toSave]

# prova 2 (più speranza ma NON VA IL MODELLO):
imageS1R1T1 <- ttt

ttt[35,35,35,1]
mask <-  test$mask

sigma <- 2
unsmooth <- GaussSmoothArray(imageS1R1T1, sigma = diag(0, 3),
                             mask = mask)
smooth <- GaussSmoothArray(imageS1R1T1,
                           sigma = diag(sigma^2, 3), mask = mask)

library(ggBrain)
plotsmooth <- ggBrain(brains = list(unsmooth, smooth),
                      mask = mask, mar = c(2,2), mar_ind = c(45,45),
                      brain_ind = c(2,2), col_ind = c("Unsmoothed", "Smoothed"),
                      type = 'signed')
plotsmooth

library(brainR)

template <- readNifti("E:/PhD/Statistical Inference II/HW5/brodmann_my.nii")

# non va piango
template1 <- niftyreg(source = brodmann, target = test,
                      scope = "affine")$image

template1[is.na(template1)] <- 0
dim(template1)

contour3d(test, level = c(1000, 2000), #add = TRUE,
          alpha = c(0.8, 0.9), color = c("yellow", "red"), mask = mask)

library(MPsychoR)
data("NeuralScanner")
onsetsS1R1 <- NeuralScanner$TRIAL_START[NeuralScanner[,2] == 1]
durationS1R1 <- NeuralScanner$RT[NeuralScanner[,2] == 1] + 1
nt <- 162
xt <- fmri.stimulus(scans = nt, onsets = onsetsS1R1,
                    dur = durationS1R1, TR = 2.5, times = TRUE)
X <- fmri.design(xt, order = 2)
head(X)

library(fields)
col5 <- colorRampPalette(c('cadetblue4', 'white', 'coral4'))
collev <- 21
maxval <- max(abs(X))
colbreak <- seq(-maxval, maxval, length.out = collev + 1)
op <- par(mar = c(5,5,5,7))
image(t(X), axes = FALSE, main = "Simple fMRI Design Matrix",
      xlab = "Contrasts", col = col5(n = collev), breaks = colbreak)
axis(1, at = c(0, 0.33, 0.66, 1),
     labels = c("E(BOLD)", "Intercept", "Linear", "Quadratic"))
box()
image.plot(t(X), legend.only = TRUE , col = col5(n = collev),
           breaks = colbreak)
par(op)

imageS1R1 <- test
dim(extractData(imageS1R1))
## [1] 53 63 46  1
spmS1R1 <- fmri.lm(imageS1R1, X)
# Error in diag(1, dy[4]) - svdresult$u %*% t(svdresult$u) : 
#   array incompatibili
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# prova 3:

# https://johnmuschelli.com/Neuroimaging_in_R/fmri_proc.html#16
library(fslr)
median_image = apply(ttt, c(1, 2, 3), median)
ortho2(median_image)

fmri.pvalue(spm, mode="basic", na.rm=FALSE, minimum.signal = 0, alpha= 0.05)

ttt


# function to calculate the t-test at each voxel and return the t value  
getT <- function(x) {    
  # can't do a t-test if variance is zero, so check before trying.  
  if (var(x) == 0) {   
    stat <- NA;   
  } else {   
    stat <- t.test(x, alternative="greater", mu=0.5)$statistic;   
  }  
  return(stat)  
}  

# plyr function aaply calls getT at each voxel in the 4d array big,  
# creating a 3d output array of t-values (what getT returns)  
t.img <- aaply(ttt, c(1,2,3), getT)









