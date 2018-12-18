

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   SPM analysis for 1 fMRI scan.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       October 2018 - November 2018
%%%%   Adapted from:  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/analysis_for_one_subject_SPM.m
%%%%                  https://github.com/wanderine/ParametricMultisubjectfMRI/blob/master/SPM/analyze_all_subjects_spm_1.m
%%%%                  https://github.com/wanderine/ParametricMultisubjectfMRI/blob/master/SPM/analyze_all_subjects_spm_checksignificance.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage      = fgetl(fopen('path_manage.txt'));
path_scratch     = fgetl(fopen('path_scratch.txt'));
HRF_models       = cellstr(['canonical   '; 'canonical_TD'; 'FIR_32_05   '; 'FIR_16_1    '; 'FIR_24_1    ']);
stim_onsets      = textscan(fopen([path_manage '/experimental_designs/stim_onsets_' subject '.txt']), '%f');
TR               = 1.97;   %-the repetition time
npts             = 261;    %-the number of time points
path_data        = [path_scratch '/scans'];

addpath('/applications/spm/spm12_7219');
addpath(genpath([path_manage '/matlab_extra_functions']));


%-initialise SPM defaults
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clear jobs;

%-data already preprocessed


for HRF_model_id = 1:length(HRF_models)

   HRF_model     = HRF_models{HRF_model_id};
   
   path_output   = [path_scratch '/analysis_output/HRF_' HRF_model];
   cd(path_output);
   
   system(['rm -r ' subject]);
   system(['mkdir ' subject]);
   
   clear jobs;
   
   %-MODEL SPECIFICATION AND ESTIMATION
   filename                                = [path_output '/' subject];
   jobs{1}.stats{1}.fmri_spec.dir          = cellstr(filename);
   jobs{1}.stats{1}.fmri_spec.timing.units = 'secs';
   jobs{1}.stats{1}.fmri_spec.timing.RT    = TR;
   
   scans = {};
   for t = 1:npts
      scans{t} = [path_data '/' subject '_SMT_wds.nii,' num2str(t)];
   end
   jobs{1}.stats{1}.fmri_spec.sess.scans        = transpose(scans);
   jobs{1}.stats{1}.fmri_spec.sess.cond(1).name = 'task1';
   
   %-specifying the experimental design
   jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = stim_onsets{1};
   jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = 0.1;
   
   %-including motion and WM/CSF covariates in the GLM
   jobs{1}.stats{1}.fmri_spec.sess.multi_reg        = {[path_manage '/additional_covariates_for_GLM/combined/add_cov_' subject '.txt']};
   
   %-specifying the HRF model, as well as the number of HRF-related covariates
   if strcmp(HRF_model, 'canonical')
      jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs   = [0 0];
      no_par                                        = 1;
   elseif strcmp(HRF_model, 'canonical_TD')
      jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs   = [1 1];
      no_par                                        = 3;
   elseif strcmp(HRF_model, 'FIR_32_05')
      jobs{1}.stats{1}.fmri_spec.bases.fir.length   = 16;
      jobs{1}.stats{1}.fmri_spec.bases.fir.order    = 32;
      no_par                                        = 32;
   elseif strcmp(HRF_model, 'FIR_16_1')
      jobs{1}.stats{1}.fmri_spec.bases.fir.length   = 16;
      jobs{1}.stats{1}.fmri_spec.bases.fir.order    = 16;
      no_par                                        = 16;
   elseif strcmp(HRF_model, 'FIR_24_1')
      jobs{1}.stats{1}.fmri_spec.bases.fir.length   = 24;
      jobs{1}.stats{1}.fmri_spec.bases.fir.order    = 24;
      no_par                                        = 24;
   end
   
   %-SPM's alternative pre-whitening: FAST
   %-following study "Accurate autocorrelation modeling substantially improves fMRI reliability":
   %-https://www.biorxiv.org/content/early/2018/09/02/323154
%   jobs{1}.stats{1}.fmri_spec.cvi              = 'FAST';

   filename_mat                                = [filename '/SPM.mat'];
   jobs{2}.stats{1}.fmri_est.spmmat            = cellstr(filename_mat);
   jobs{3}.stats{1}.con.spmmat                 = cellstr(filename_mat);
   jobs{3}.stats{1}.con.consess{1}.fcon        = struct('name', 'task1 > rest', 'convec', eye(no_par), 'sessrep', 'none');
   jobs{4}.stats{1}.results.spmmat             = cellstr(filename_mat);
   jobs{4}.stats{1}.results.conspec.contrasts  = 1;
   jobs{4}.stats{1}.results.conspec.threshdesc = 'none';
   jobs{4}.stats{1}.results.conspec.thresh     = 0.001;
   jobs{4}.stats{1}.results.conspec.extent     = 0;
   jobs{4}.stats{1}.results.print              = false;

   spm_jobman('run', jobs);
   clear jobs;

   %-MULTIPLE COMPARISON CORRECTION

   cd(filename);
   
   %-get F-map
   V              = spm_vol([path_output '/' subject '/spmF_0001.nii']);
   [tmap, aa]     = spm_read_vols(V);

   %-calculate cluster extent threshold
   %-using the default cluster defining threshold (CDT) of 0.001
   [k, Pc]        = corrclusth(SPM, 0.001, 0.05, 1:100000);

   %-get the degrees of freedom ('df1' and 'df2')
   df             = [no_par SPM.xX.erdf];

   u              = spm_u(0.001, df, 'F');
   indices        = find(tmap > u);
   fid            = fopen('indices.txt', 'wt');
   fprintf(fid, '%8.0f\n', indices);
   
   %-get the size of the largest cluster
   max_cluster    = max_extent(tmap, indices);
   if max_cluster >= k
      indices     = find(tmap>u);
   else
      indices     = '';
   end
   fid            = fopen('indices.txt', 'wt');
   fprintf(fid, '%8.0f\n', indices);
   
   clear jobs;
   jobs{1}.stats{1}.results.spmmat             = cellstr(filename_mat);
   jobs{1}.stats{1}.results.conspec.contrasts  = 1;
   jobs{1}.stats{1}.results.conspec.threshdesc = 'none';
   jobs{1}.stats{1}.results.conspec.thresh     = 0.001;
   jobs{1}.stats{1}.results.conspec.extent     = 0;
   jobs{1}.stats{1}.results.print              = true;

   spm_jobman('run', jobs);

   %-saving 3D maps showing where the significant clusters are located
   cd([path_output '/' subject]);
   mask_MNI                    = niftiread('mask.nii');
   cluster_binary_MNI          = zeros(size(mask_MNI));
   indices                     = fscanf(fopen('indices.txt', 'r'), '%d');
   cluster_binary_MNI(indices) = 1;
   niftiwrite(cluster_binary_MNI, 'cluster_binary_MNI.nii');
   
end

cd(path_manage)
