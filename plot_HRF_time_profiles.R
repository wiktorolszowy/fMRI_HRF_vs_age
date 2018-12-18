

###############################################################################################
####   Analysing the impact of age on HRF.
####   Written by:    Wiktor Olszowy, University of Cambridge
####   Contact:       wo222@cam.ac.uk
####   Created:       July 2018 - December 2018
###############################################################################################


library(AnalyzeFMRI)
library(R.matlab)
library(parallel)
library(plotrix)    #-for color legend bars

path_manage   = readLines("path_manage.txt")
path_scratch  = readLines("path_scratch.txt")
subjects      = list.files(paste0(path_scratch, "/scans"))
subjects      = unique(substr(subjects, 1, 8))
setwd(path_manage)
data          = read.csv("CamCAN_subjects_info.csv", sep=",")
anat          = f.read.nifti.volume("brain_parcellation/rlabels_Neuromorphometrics.nii")
#-does not work with ||, needs to be |
A1_loc        = which(anat==200 | anat==201)
M1_loc        = which(anat==182 | anat==183 | anat==192 | anat==193)
V1_loc        = which(anat==108 | anat==109)
VOIs          = c(200, 201, 108, 109, 182, 183, 192, 193)
VOIs_t        = c("Right superior temporal gyrus", "Left superior temporal gyrus", "Right calcarine cortex", "Left calcarine cortex", "Right precentral gyrus", "Left precentral gyrus", "Right supplementary motor cortex", "Left supplementary motor cortex")
combined_perc = readMat("combined_results/combined_perc.mat")
combined_perc = combined_perc[[1]]

combined_perc[which(combined_perc==-1)] = NA

#-altogether
no_subjects = length(subjects)

age     = rep(-1, no_subjects)
gender  = rep(data$Gender[1], no_subjects)
for (i in 1:no_subjects) {
   age[i]    = data$Age   [which(paste0("sub-", subjects[i])==paste0("sub-", data$CCID))]
   gender[i] = data$Gender[which(paste0("sub-", subjects[i])==paste0("sub-", data$CCID))]
}

#-sorting subjects w.r.t. age (youngest to oldest)
subjects = subjects[order(age)]

path_output_top  = paste0(path_scratch, "/analysis_output")

setwd(path_manage)


HRF_canonical_DD = read.table(paste0(path_manage, "/HRF_canonical_TD/HRF_canonical_TD.csv"), header=F)
canonical        = HRF_canonical_DD$V1
temp_der         = HRF_canonical_DD$V2
disp_der         = HRF_canonical_DD$V3
time_scale       = seq(0, 32, by=0.1)


#-save subjects for which no <0.001 activation in A1 or V1
A1_started = F
V1_started = F
for (subject_id in 1:no_subjects) {
   HRF_model      = "canonical_TD"
   subject        = subjects[subject_id]
   path           = paste0(path_output_top, "/HRF_", HRF_model, "/", subject)
   if (file.exists(paste0(path, "/spmF_0001.nii"))) {
      setwd(path)
      F_stat      = f.read.nifti.volume("spmF_0001.nii")
      #-dfs taken from header, the same in SPM.xX.erdf
      thr_F       = qf(0.999, 3, 243)
      setwd(path_manage)
      if (sum(F_stat[which(anat==200 | anat==201)] > thr_F)==0) {
         if (A1_started==F) {
            write(substr(subject, nchar(subject)-5, nchar(subject)), file="no_activation_in_A1.txt")
         } else {
            write(substr(subject, nchar(subject)-5, nchar(subject)), file="no_activation_in_A1.txt", append=T)
         }
         A1_started = T
      }
      if (sum(F_stat[which(anat==108 | anat==109)] > thr_F)==0) {
         if (V1_started==F) {
            write(substr(subject, nchar(subject)-5, nchar(subject)), file="no_activation_in_V1.txt")
         } else {
            write(substr(subject, nchar(subject)-5, nchar(subject)), file="no_activation_in_V1.txt", append=T)
         }
         V1_started = T
      }
   }
}


#-calculating average HRF for 'canonical_TD'
HRF_est_canTD_list = mclapply(1:no_subjects, function(subject_id) {
   HRF_model      = "canonical_TD"
   subject        = subjects[subject_id]
   cat(subject, "\n")
   HRF_est        = array(NA, dim=c(length(VOIs), length(time_scale)))
   path           = paste0(path_output_top, "/HRF_", HRF_model, "/", subject)
   if (file.exists(paste0(path, "/spmF_0001.nii"))) {
      setwd(path)
      beta_1      = f.read.nifti.volume("beta_0001.nii")
      beta_2      = f.read.nifti.volume("beta_0002.nii")
      beta_3      = f.read.nifti.volume("beta_0003.nii")
      F_stat      = f.read.nifti.volume("spmF_0001.nii")
      #-dfs taken from header, the same in SPM.xX.erdf
      thr_F       = qf(0.999, 3, 243)
      for (i in 1:length(VOIs)) {
         VOI         = VOIs[i]
         beta_1_av   = mean(beta_1[which(F_stat>thr_F & anat==VOI & !is.na(beta_1) & !is.na(beta_2) & !is.na(beta_3))])
         beta_2_av   = mean(beta_2[which(F_stat>thr_F & anat==VOI & !is.na(beta_1) & !is.na(beta_2) & !is.na(beta_3))])
         beta_3_av   = mean(beta_3[which(F_stat>thr_F & anat==VOI & !is.na(beta_1) & !is.na(beta_2) & !is.na(beta_3))])
         HRF_est[i,] = beta_1_av %*% t(canonical) + beta_2_av %*% t(temp_der) + beta_3_av %*% t(disp_der)
      }
   }
   return(HRF_est)
}, mc.cores=24)
HRF_est_canTD = array(NA, dim=c(no_subjects, length(VOIs), length(time_scale)))
for (subject_id in 1:no_subjects) {
   HRF_est_canTD[subject_id, , ] = HRF_est_canTD_list[[subject_id]]
}
HRF_est_canTD_across_sub = apply(HRF_est_canTD, c(2,3), mean, na.rm=T)


#-calculating average HRF for 'FIR_32_05'
params_FIR_list = mclapply(1:no_subjects, function(subject_id) {
   HRF_model       = "FIR_32_05"
   subject         = subjects[subject_id]
   cat(subject, "\n")
   betas_means     = array(NA, c(32, length(VOIs)))
   path            = paste0(path_output_top, "/HRF_", HRF_model, "/", subject)
   if (file.exists(paste0(path, "/spmF_0001.nii"))) {
      setwd(path)
      beta            = array(NA, dim=c(32, length(as.vector(f.read.nifti.volume("beta_0001.nii")))))
      leave_ids_prev  = 1:length(as.vector(f.read.nifti.volume("beta_0001.nii")))
      F_stat          = f.read.nifti.volume("spmF_0001.nii")
      #-dfs taken from header, the same in SPM.xX.erdf
      thr_F           = qf(0.999, 32, 214)
      for (par_id in 1:32) {
         if (par_id < 10) {
            beta[par_id,] = as.vector(f.read.nifti.volume(paste0("beta_000", par_id, ".nii")))
         } else {
            beta[par_id,] = as.vector(f.read.nifti.volume(paste0("beta_00",  par_id, ".nii")))
         }
         leave_ids        = which(!is.na(beta[par_id,]) & beta[par_id,]!=0)
         leave_ids        = intersect(leave_ids, leave_ids_prev)
         leave_ids_prev   = leave_ids
      }
      for (i in 1:length(VOIs)) {
         VOI = VOIs[i]
         voxels_set       = intersect(intersect(which(F_stat>thr_F), which(anat==VOI)), leave_ids)
         if (length(voxels_set)==1) {
            betas_means[, i] = beta[, voxels_set]
         } else {
            betas_means[, i] = apply(beta[, voxels_set], 1, mean)
         }
      }
   }
   return(betas_means)
}, mc.cores=24)
params_FIR = array(NA, dim=c(no_subjects, 32, length(VOIs)))
for (subject_id in 1:no_subjects) {
   params_FIR[subject_id, , ] = params_FIR_list[[subject_id]]
}
#-TR of 1.97s, 32 bins covering 32*0.5s=16s, 1 regressor is assigned to each of the 32 bins
#-FIR is not a continuous function (no interpolation opposed to 'tent' in AFNI), thus /2!
time_FIR      = seq(0.5/2, 0.5/2+(32-1)*0.5, by=0.5)
response_FIR  = apply(params_FIR, c(2,3), mean, na.rm=T)
HRF_est_FIR   = array(NA, dim=dim(HRF_est_canTD))
for (subject_id in 1:no_subjects) {
   for (i in 1:length(VOIs)) {
      for (time_bin in 1:dim(params_FIR)[2]) {
         HRF_est_FIR[subject_id, i, ((time_bin-1)*5+1) : (time_bin*5)] = rep(params_FIR[subject_id, time_bin, i], 5)
      }
   }
}
HRF_est_FIR[which(is.nan(HRF_est_FIR))] = NA


#-calculating % of active voxels to all voxels within each VOI for each subject
perc_act_list   = mclapply(1:no_subjects, function(subject_id) {
   subject      = subjects[subject_id]
   perc_act_aux = array(NA, dim=c(2, length(VOIs)))
   for (HRF_model_id in 1:2) {
      if (HRF_model_id == 1) {
         HRF_model = "canonical_TD"
         thr_F     = qf(0.999, 3,  243)
      } else {
         HRF_model = "FIR_32_05"
         thr_F     = qf(0.999, 32, 214)
      }
      setwd(paste0(path_output_top, "/HRF_", HRF_model, "/", subject))
      if (file.exists("spmF_0001.nii")) {
         F_stat  = f.read.nifti.volume("spmF_0001.nii")
         #-dfs taken from header, the same in SPM.xX.erdf
         for (i in 1:length(VOIs)) {
            VOI  = VOIs[i]
            #-in percentages
            #-as the VOIs non-empty, there is no division by zero
            perc_act_aux[HRF_model_id, i] = 100*length(intersect(which(F_stat>thr_F), which(anat==VOI))) / length(which(anat==VOI))
         }
      } else {
         cat(getwd(), "\n")
      }
   }
   return(perc_act_aux)
}, mc.cores=24)
perc_act = array(NA, dim=c(no_subjects, 2, length(VOIs)))
for (subject_id in 1:no_subjects) {
   perc_act[subject_id, , ] = perc_act_list[[subject_id]]
}


setwd(path_manage)
pdf("figures/HRF_age_time_profiles_across_VOIs.pdf", width=8, height=11)

   #-horizontal y axis labels
   par(las=1)

   par(mfrow=c(2,1))
   matplot(time_scale, t(HRF_est_canTD_across_sub), type="l", col=rainbow(8), lwd=2, lty=1, xlim=c(0,16), ylim=c(-0.15, 0.2), xlab="Post-stimulus time [s]", ylab="", main="Estimated HRF for canonical + TD")
   abline(v=5, col="gray")
   abline(h=0, col="gray")
   legend("topright", legend=VOIs_t, col=rainbow(8), ncol=1, lty=1, lwd=4.5, bty="n")
   matplot(time_FIR, response_FIR, type="l", col=rainbow(8), lwd=2, lty=1, xlim=c(0,16), ylim=c(-0.15, 0.2), xlab="Post-stimulus time [s]", ylab="", main="Estimated HRF for FIR (32 x 0.5s bins)")
   abline(v=5, col="gray")
   abline(h=0, col="gray")
   legend("topright", legend=VOIs_t, col=rainbow(8), ncol=1, lty=1, lwd=4.5, bty="n")

dev.off()


#-from https://stackoverflow.com/questions/14660372/common-main-title-of-a-figure-panel-compiled-with-parmfrow
line2user = function(line, side) {
   lh    = par('cin')[2] * par('cex') * par('lheight')
   x_off = diff(grconvertX(0:1, 'inches', 'user'))
   y_off = diff(grconvertY(0:1, 'inches', 'user'))
   switch(side,
      `1` = par('usr')[3] - line * y_off * lh,
      `2` = par('usr')[1] - line * x_off * lh,
      `3` = par('usr')[4] + line * y_off * lh,
      `4` = par('usr')[2] + line * x_off * lh,
      stop("side must be 1, 2, 3, or 4", call.=FALSE))
}


#-from https://stackoverflow.com/questions/25015410/r-find-full-width-at-half-maximum-for-a-gausian-density-distribution
fwhm_fun = function(vec) {
   if (length(na.omit(vec)) < 2) {
      return (NA)
   }
   d     = density(na.omit(vec))
   xmax  = d
   xmax  = d$x[d$y==max(d$y)]
   x1    = d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax]-max(d$y)/2))]
   x2    = d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax]-max(d$y)/2))]
   return (x2-x1)
}


pdf("figures/HRF_age_time_profiles_shape_features.pdf", width=6, height=10)

   cex = 0.4
   lwd = 2

   #-horizontal y axis labels
   par(las=1)

   #-fills the matrix by columns
   par(mfcol=c(6,2))

   ########################### canTD  #################################################################

   #-remove outliers
   for (i in 1:length(VOIs)) {
      for (subject_id in 1:no_subjects) {
         where_maximum = which.max(HRF_est_canTD[subject_id, i, ])
         if (length(where_maximum)>0) {
            peak_time  = time_scale[where_maximum]
            if (peak_time < 2 || peak_time > 8) {
               perc_act[subject_id, 1, i]     = NA
               HRF_est_canTD[subject_id, i, ] = NA
            }
         }
      }
   }

   #-just to create space for figure titles
   plot.new()

   par(mgp = c(1.75, 0.5, 0))
   #-c(bottom, left, top, right)
   par(mai = c(0.2, 0.4, 0.25, 0.05))

   #-perc of sig voxels
   matplot(age, perc_act[, 1, ], pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(perc_act[, 1, ], na.rm=T)))
   for (i in 1:length(VOIs)) {
      out = which(is.na(perc_act[, 1, i]))
      if (length(out) > 0) {
         points(age[-out], lm(perc_act[, 1, i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
      } else {
         points(age,       lm(perc_act[, 1, i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
      }
   }
   mtext("canonical + TD", font=2, side=3, line=2, cex=0.8)

   #-maximum
   maximum = apply(HRF_est_canTD, c(1,2), max, na.rm=F)
   matplot(age, maximum, pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(maximum, na.rm=T)))
   for (i in 1:length(VOIs)) {
      out = which(is.na(maximum[,i]))
      points(age[-out], lm(maximum[,i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
   }

   #-peak time
   peak_time = array(NA, dim=dim(maximum))
   for (i in 1:length(VOIs)) {
      for (subject_id in 1:no_subjects) {
         where_maximum = which.max(HRF_est_canTD[subject_id, i, ])
         if (length(where_maximum)>0) {
            peak_time[subject_id, i] = time_scale[where_maximum]
         }
      }
   }
   matplot(age, peak_time, pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(peak_time, na.rm=T)))
   abline(h=5, lwd=lwd, col="darkgray")
   for (i in 1:length(VOIs)) {
      out      = which(is.na(maximum[,i]))
      outliers = which(peak_time < 2 || peak_time > 9)
      points(age[-out], lm(peak_time[,i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
   }

   #-width
   width = apply(HRF_est_canTD, c(1,2), fwhm_fun)
   matplot(age, width, pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(width, na.rm=T)))
   abline(h=5, lwd=lwd, col="darkgray")
   for (i in 1:length(VOIs)) {
      out = which(is.na(maximum[,i]))
      points(age[-out], lm(width[,i]~age)$fitted,     col=rainbow(8)[i], type="l", lwd=lwd)
   }

   if (i==4 || i==8) {
      mtext("Post-stimulus time [s]", side=1, line=2, cex=0.75)
   }

   plot.new()
   legend(0.22, 0.98, legend=VOIs_t[1:4], col=rainbow(8)[1:4], ncol=1, lty=1, lwd=lwd, bty="n")

   ########################### FIR  ###################################################################

   #-remove outliers
   for (i in 1:length(VOIs)) {
      for (subject_id in 1:no_subjects) {
         where_maximum = which.max(HRF_est_FIR[subject_id, i, ])
         if (length(where_maximum)>0) {
            peak_time  = time_scale[where_maximum]
            if (peak_time < 2 || peak_time > 8) {
               perc_act[subject_id, 2, i]   = NA
               HRF_est_FIR[subject_id, i, ] = NA
            }
         }
      }
   }

   maximum   = array(NA, dim=c(no_subjects, length(VOIs)))
   peak_time = array(NA, dim=c(no_subjects, length(VOIs)))
   width     = array(NA, dim=c(no_subjects, length(VOIs)))
   for (subject_id in 1:no_subjects) {
      for (i in 1:length(VOIs)) {
         if (sum(is.na(HRF_est_FIR[subject_id, i, ])) < length(time_scale)) {
            maximum_single = max(HRF_est_FIR[subject_id, i, ], na.rm=T)
            time_point     = median(which(HRF_est_FIR[subject_id, i, ] == maximum_single))
            maximum[subject_id, i]   = maximum_single
            peak_time[subject_id, i] = time_scale[time_point]
            width[subject_id, i]     = fwhm_fun(HRF_est_FIR[subject_id, i, ])
         } else {
            maximum[subject_id, i]   = NA
            peak_time[subject_id, i] = NA
            width[subject_id, i]     = NA
         }
      }
   }

   #-just to create space for figure titles
   plot.new()

   #-perc of sig voxels
   matplot(age, perc_act[, 2, ], pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(perc_act[, 2, ], na.rm=T)))
   for (i in 1:length(VOIs)) {
      out = which(is.na(perc_act[, 2, i]))
      if (length(out) > 0) {
         points(age[-out], lm(perc_act[, 2, i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
      } else {
         points(age,       lm(perc_act[, 2, i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
      }
   }
   mtext("FIR (32 x 0.5s bins)", font=2, side=3, line=2, cex=0.8)

   text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
        line2user(line=1, side=3), expression(bold("% of p<0.001 voxels within the VOI mask")), xpd=NA)

   #-maximum
   matplot(age, maximum, pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(maximum, na.rm=T)))
   for (i in 1:length(VOIs)) {
      out = which(is.na(maximum[,i]))
      points(age[-out], lm(maximum[,i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
   }
   text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
        line2user(line=1, side=3), expression(bold("HRF maximum [a.u.]")), xpd=NA)

   #-peak time
   matplot(age, peak_time, pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(peak_time, na.rm=T)))
   abline(h=5, lwd=lwd, col="darkgray")
   for (i in 1:length(VOIs)) {
      out = which(is.na(maximum[,i]))
      points(age[-out], lm(peak_time[,i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
   }
   text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
        line2user(line=1, side=3), expression(bold("HRF peak time [s]")), xpd=NA)

   #-width
   matplot(age, width, pch=16, cex=cex, xlab="subject's age", ylab="", col=rainbow(8), ylim=c(0, max(width, na.rm=T)))
   abline(h=5, lwd=lwd, col="darkgray")
   for (i in 1:length(VOIs)) {
      out = which(is.na(maximum[,i]))
      points(age[-out], lm(width[,i]~age)$fitted, col=rainbow(8)[i], type="l", lwd=lwd)
   }
   text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
        line2user(line=1, side=3), expression(bold("HRF width [s]")), xpd=NA)

   if (i==4 || i==8) {
      mtext("Post-stimulus time [s]", side=1, line=2, cex=0.75)
   }

   plot.new()
   legend(0, 0.98, legend=VOIs_t[5:8], col=rainbow(8)[5:8], ncol=1, lty=1, lwd=lwd, bty="n")

dev.off()


for (smooth in c(F, T)) {

   for (i in 1:length(VOIs)) {

      if (i==1) {
         if (smooth==T) {
            pdf("figures/HRF_age_time_profiles_smoothed_1.pdf",   width=6, height=10.7)
         } else {
            pdf("figures/HRF_age_time_profiles_unsmoothed_1.pdf", width=6, height=10.7)
         }
         #-horizontal y axis labels
         par(las=1)
         #-'mfrow=c(6, 2)' used as a trick so that one blank subplot above and one blank subplot below available
         #-this is necessary as the margins are very small, and for the first and last real subplots there are elements beyond margins
         par(mfrow=c(6, 2))
         #-just to create space for figure titles
         plot.new()
         plot.new()
      } else if (i==5) {
         dev.off()
         if (smooth==T) {
            pdf("figures/HRF_age_time_profiles_smoothed_2.pdf",   width=6, height=10.7)
         } else {
            pdf("figures/HRF_age_time_profiles_unsmoothed_2.pdf", width=6, height=10.7)
         }
         #-horizontal y axis labels
         par(las=1)
         #-'mfrow=c(6, 2)' used as a trick so that one blank subplot above and one blank subplot below available
         #-this is necessary as the margins are very small, and for the first and last real subplots there are elements beyond margins
         par(mfrow=c(6, 2))
         #-just to create space for figure titles
         plot.new()
         plot.new()
      }

      par(mgp = c(1.75, 0.5, 0))
      par(mai = c(0.2, 0.4, 0.23, 0.05))

      map_canTD   = HRF_est_canTD[, i, ]
      map_FIR     = HRF_est_FIR  [, i, ]
      map_FIR[,1] = 0
      out_canTD   = c()
      out_FIR     = c()

      for (subject_id in 1:no_subjects) {
         if (sum(is.na(map_canTD[subject_id, ])) > 0) {
            out_canTD = c(out_canTD, subject_id)
         }
         if (is.na(map_FIR[subject_id, 2])) {
            out_FIR   = c(out_FIR, subject_id)
         } else {
            map_FIR[subject_id, is.na(map_FIR[subject_id, ])] = 0
         }
      }

      #-canTD
      z_val = t(map_canTD[-out_canTD, ])
      if (smooth==T) {
         for (time_point in 1:dim(z_val)[1]) {
            #-bandwidth: 5 subjects
            z_val[time_point, ] = ksmooth(1:length(z_val[time_point, ]), z_val[time_point, ], n.points=length(z_val[time_point, ]), kernel="normal", bandwidth=5)$y
         }
      }
      #-truncating the z matrix, so that no white color in the heatmaps (below and above zlim)
      z_lims    = 1.2*c(-0.07, 0.09)
      z_val_aux = z_val
      z_val_aux[which(z_val < z_lims[1])] = z_lims[1]
      z_val_aux[which(z_val > z_lims[2])] = z_lims[2]
      image(x=1:321, y=1:length(age[-out_canTD]), z=z_val_aux, col=rev(rainbow(20)), xlab="Post-stimulus time [s]", ylab="Age [years]", xaxt="n", yaxt="n", xlim=c(0,170), zlim=z_lims)
      age_in = age[-out_canTD]
      axis(1, at=seq(1,  301,by=50),  labels=time_scale[seq(1,301,by=50)])
      axis(2, at=seq(100,600,by=100), labels=age_in[seq(100,600,by=100)])
      abline(v=seq(51,301,by=50), lty="dashed")
      if (i==1 || i==5) {
         mtext("canonical + TD", font=2, side=3, line=2, cex=0.8)
      }
      if (i==4 || i==8) {
         mtext("Post-stimulus time [s]", side=1, line=2, cex=0.75)
      }
      
      #-FIR
      z_val = t(map_FIR[-out_FIR, ])
      if (smooth==T) {
         for (time_point in 1:dim(z_val)[1]) {
            #-bandwidth: 5 subjects
            z_val[time_point, ] = ksmooth(1:length(z_val[time_point, ]), z_val[time_point, ],  n.points=length(z_val[time_point, ]), kernel="normal", bandwidth=5)$y
         }
      }
      #-additional smoothing for OHBM abstract (only FIR!!!!)
#      if (smooth==T) {
#         for (subject_id in 1:dim(z_val)[2]) {
#            #-bandwidth: 7*0.1s = 0.7s
#            z_val[, subject_id] = ksmooth(1:length(z_val[, subject_id]), z_val[, subject_id ], n.points=length(z_val[, subject_id]), kernel="normal", bandwidth=7)$y
#         }
#      }
      #-truncating the z matrix, so that no white color in the heatmaps (below and above zlim)
      z_lims    = 1.2*c(-0.12, 0.09/0.07*0.12)
      z_val_aux = z_val
      z_val_aux[which(z_val < z_lims[1])] = z_lims[1]
      z_val_aux[which(z_val > z_lims[2])] = z_lims[2]
      image(x=1:321, y=1:length(age[-out_FIR]),   z=z_val_aux, col=rev(rainbow(20)), xlab="Post-stimulus time [s]", ylab="Age [years]", xaxt="n", yaxt="n", xlim=c(0,170), zlim=z_lims)
      age_in = age[-out_FIR]
      axis(1, at=seq(1,  301,by=50),  labels=time_scale[seq(1,301,by=50)])
      axis(2, at=seq(100,600,by=100), labels=age_in[seq(100,600,by=100)])
      abline(v=seq(51,301,by=50), lty="dashed")
      if (i==1 || i==5) {
         mtext("FIR (32 x 0.5s bins)", font=2, side=3, line=2, cex=0.8)
      }
      if (i==4 || i==8) {
         mtext("Post-stimulus time [s]",       side=1, line=2, cex=0.75)
      }

      text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
           line2user(line=1, side=3), bquote(bold(.(VOIs_t[i]))), xpd=NA)

      #-adding color legend bar
      if (i==4 || i==8) {
         #-canTD
         plot.new()
         z_lims = 1.2*c(-0.07, 0.09)
         color.legend(0.1, 0.81, 0.9, 0.9, legend=c(round(z_lims[1], 2), rep("", 8), 0, rep("", 9), round(z_lims[2], 2)), rect.col=rev(rainbow(20)), align="lt", gradient="x", cex=0.75)
         #-FIR
         plot.new()
         z_lims = 1.2*c(-0.12, 0.09/0.07*0.12)
         color.legend(0.1, 0.81, 0.9, 0.9, legend=c(round(z_lims[1], 2), rep("", 8), 0, rep("", 9), round(z_lims[2], 2)), rect.col=rev(rainbow(20)), align="lt", gradient="x", cex=0.75)
      }

   }

   dev.off()

}


#-removing space around figures
setwd(paste0(path_manage, "/figures"))
system("pdfcrop HRF_age_time_profiles_across_VOIs.pdf    HRF_age_time_profiles_across_VOIs.pdf")
system("pdfcrop HRF_age_time_profiles_shape_features.pdf HRF_age_time_profiles_shape_features.pdf")
system("pdfcrop HRF_age_time_profiles_smoothed_1.pdf     HRF_age_time_profiles_smoothed_1.pdf")
system("pdfcrop HRF_age_time_profiles_smoothed_2.pdf     HRF_age_time_profiles_smoothed_2.pdf")
system("pdfcrop HRF_age_time_profiles_unsmoothed_1.pdf   HRF_age_time_profiles_unsmoothed_1.pdf")
system("pdfcrop HRF_age_time_profiles_unsmoothed_2.pdf   HRF_age_time_profiles_unsmoothed_2.pdf")
setwd(path_manage)
