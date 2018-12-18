

#-Wiktor Olszowy

library(R.matlab)

path_manage  = readLines("path_manage.txt")
path_scratch = readLines("path_scratch.txt")
npts         = 261
TR           = 1.97
HRF_models   = c("canonical", "canonical_TD", "FIR_32_05", "FIR_16_1", "FIR_24_1")


#-copy pre-processed scans
library(parallel)
path_to_preprocessed = "/home/wo222/scratch/fMRI_data/CamCAN_preprocessed/aamod_waveletdespike_00001"
path_destination     = "/home/wo222/scratch/fMRI_method_validations/modelling_age_HRF/scans"
setwd(path_to_preprocessed)
subjects             = list.files()
subjects_out         = c("CC610462", "CC710518", "CC110045", "CC120065", "CC221755", "CC410129", "CC410222", "CC510639", "CC610050", "CC610146", "CC710214")
subjects             = subjects[-which(subjects %in% subjects_out)]
#res = mclapply(1:length(subjects), function(subject_id) {
#   subject = subjects[subject_id]
#   setwd(paste0(path_to_preprocessed, "/", subject, "/SMT/"))
#   #-copy the wavelet-despiked scan
#   wds = list.files()[which(substr(list.files(), nchar(list.files())-7, nchar(list.files()))=="_wds.nii")]
#   system(paste0("cp ", path_to_preprocessed, "/", subject, "/SMT/", wds, " ", path_destination, "/", subject, "_SMT_wds.nii"))
#   #-copy the stimulus timing file
#   system(paste0("cp /home/wo222/scratch/fMRI_data/CamCAN/cc700/mri/pipeline/release004/BIDSsep/func_smt/sub-", subject, "/func/sub-", subject, "_task-SMT_events.tsv ", path_destination, "/", subject, "_SMT_events.tsv"))
#}, mc.cores=24)

#-make folders
setwd(path_scratch)
system("mkdir analysis_output")
setwd(paste0(path_scratch, "/analysis_output"))
for (HRF_model_id in 1:length(HRF_models)) {
   system(paste0("mkdir HRF_", HRF_models[HRF_model_id]));
}

#-make experimental designs
setwd(path_manage)
system("mkdir experimental_designs")
setwd(paste0(path_manage, "/experimental_designs"))
for (subject in subjects) {
   setwd(paste0(path_scratch, "/scans"))
   events_times = read.table(file=paste0(subject, "_SMT_events.tsv"), sep="\t", header=T)
   onsets       = events_times$onset
   setwd(paste0(path_manage, "/experimental_designs"))
   write(events_times$onset, file=paste0("stim_onsets_", subject, ".txt"), sep="\n")
}

#-make parallel commands
setwd(path_manage)
for (subject_id in 1:length(subjects)) {
   subject = subjects[subject_id]
   cat("matlab -r -nodesktop \"subject='", subject, "'; run('analysis_for_one_subject_SPM.m'); exit\" \n", file=paste0("parallel_commands/command_1_", subject_id, ".sh"), sep="", append=F)
}
for (subject_id in 1:length(subjects)) {
   subject = subjects[subject_id]
   cat("matlab -r -nodesktop \"subject='", subject, "'; run('make_VOIs.m'); exit\" \n",                    file=paste0("parallel_commands/command_2_", subject_id, ".sh"), sep="", append=F)
}
i = 0;
for (HRF_model_id in 1:6) {
   for (age in 1:2) {
      i = i + 1;
      cat("matlab -r -nodesktop \"HRF_model_id=", HRF_model_id, "; age=", age, ";  run('make_group_analyses.m'); exit\" \n", file=paste0("parallel_commands/command_3_", i, ".sh"), sep="", append=F)
   }
}

subjects_all = read.csv("CamCAN_subjects_info.csv")
write.table(subjects_all$CCID, "subjects_CCID_all_all.txt", row.names=F, col.names=F, quote=F)

#-for each subject, combine motion covariates and WM/CSF covariates to one '.txt' file
system("mkdir additional_covariates_for_GLM/combined")
for (subject in list.files(paste0(path_manage, "/additional_covariates_for_GLM/motion/"))) {
   motion_filename = list.files(paste0(path_manage, "/additional_covariates_for_GLM/motion/", subject, "/SMT/"))
   motion          = as.matrix(read.table(paste0(path_manage, "/additional_covariates_for_GLM/motion/", subject, "/SMT/", motion_filename)))
   #-for one subject, CBU data did not include WM/CSF covariates
   if (subject != 'CC721052') {
      WM_and_CSF      = readMat(paste0(path_manage, "/additional_covariates_for_GLM/WM_and_CSF/", subject, "/SMT/compSignal.mat"))
      combined        = matrix(NA, nrow=261, ncol=8)
      combined[,1:6]  = motion;
      combined[,7:8]  = WM_and_CSF$compTC[,2:3];
   } else {
      combined        = matrix(NA, nrow=261, ncol=6)
      combined[,1:6]  = motion;
   }
   write.table(format(combined, digits=8), paste0(path_manage, "/additional_covariates_for_GLM/combined/add_cov_", subject, ".txt"), col.names=F, row.names=F, sep="\t", quote=F)
}

setwd(path_manage)
