

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Creating VOIs. Needed for the later balloon model estimation.
%%%%   Written by:    Rik Henson + Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       July 2018 - December 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage  = fgetl(fopen('path_manage.txt'));
path_scratch = fgetl(fopen('path_scratch.txt'));
HRF_models   = cellstr(['canonical   '; 'canonical_TD'; 'FIR_32_05   '; 'FIR_16_1    '; 'FIR_24_1    ']);
VOIs         = cellstr(['200'; '201'; '108'; '109'; '182'; '183'; '192'; '193']);

%-VOI numbers refer to /applications/spm/spm12_7219/tpm/labels_Neuromorphometrics.nii


addpath('/applications/spm/spm12_7219');

%-initialise SPM defaults
spm('Defaults', 'fMRI');
spm_jobman('initcfg');


for HRF_model_id = 1:length(HRF_models)

   HRF_model     = HRF_models{HRF_model_id};

   %-specifying the number of HRF-related covariates
   if strcmp(HRF_model, 'canonical')
      no_par = 1;
   elseif strcmp(HRF_model, 'canonical_TD')
      no_par = 3;
   elseif strcmp(HRF_model, 'FIR_32_05')
      no_par = 32;
   elseif strcmp(HRF_model, 'FIR_16_1')
      no_par = 16;
   elseif strcmp(HRF_model, 'FIR_24_1')
      no_par = 24;
   end

   for i  = 1:length(VOIs)

      VOI = VOIs{i};

      %-different HRF models only for left STG
      if ~strcmp(VOI, '201') && ~strcmp(HRF_model, 'canonical_TD')
         continue;
      end

      path_output_top = [path_scratch '/analysis_output/HRF_' HRF_model '/' subject];

      if exist([path_output_top '/cluster_binary_MNI.nii'], 'file')

         cd(path_output_top);
         system(['rm -rf VOI_' VOI]);
         system(['mkdir  VOI_' VOI]);
         path_output = [path_output_top '/VOI_' VOI];
         system(['cp *.* VOI_' VOI '/']);
         cd(path_output);

         load SPM;

         SPM.swd  = path_output;
         SPM.xCon = spm_FcUtil('Set', 'EOI', 'F', 'c', [eye(no_par) zeros(no_par, size(SPM.xX.X,2)-no_par)]', SPM.xX.xKXs);

         %-why do I need it?
         spm_contrasts(SPM, 1);

         base  = [path_scratch '/scans/' subject '_SMT_wds.nii'];
         tmp   = {};
         for f = 1:length(SPM.xY.P)
            tmp{f} = [base sprintf(',%d', f)];
         end

         SPM.xY.P  = strvcat(tmp);
         SPM.xY.VY = spm_vol(SPM.xY.P);

         save SPM SPM;

         %-ROI raw data

         clear matlabbatch;

         SPM_file = [path_output '/SPM.mat'];

         matlabbatch{1}.spm.util.voi.spmmat  = cellstr(SPM_file);
         matlabbatch{1}.spm.util.voi.adjust  = 1;
         matlabbatch{1}.spm.util.voi.session = 1;
         matlabbatch{1}.spm.util.voi.name    = VOI;

         str = [path_manage '/brain_parcellation/mask_VOI_' VOI '.nii,1'];

         matlabbatch{1}.spm.util.voi.roi{1}.mask.image      = cellstr(str);
         matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold  = 0.5;
         matlabbatch{1}.spm.util.voi.expression             = 'i1';

         matlabbatch{1}.spm.util.voi.roi{2}.spm.spmmat      = {''}; %cellstr(SPM_file);
         matlabbatch{1}.spm.util.voi.roi{2}.spm.contrast    = 1;
         matlabbatch{1}.spm.util.voi.roi{2}.spm.conjunction = 1;
         matlabbatch{1}.spm.util.voi.roi{2}.spm.threshdesc  = 'none';
         matlabbatch{1}.spm.util.voi.roi{2}.spm.thresh      = 0.001;
         matlabbatch{1}.spm.util.voi.roi{2}.spm.extent      = 0;
         matlabbatch{1}.spm.util.voi.roi{2}.spm.mask        = struct('contrast', {}, 'thresh', {}, 'mtype', {});
         matlabbatch{1}.spm.util.voi.expression             = 'i1&i2';

         spm_jobman('run', matlabbatch);

         tmp = SPM;
         clear SPM;

         SPM.Sess  = tmp.Sess;
         SPM.xY.RT = tmp.xY.RT;

         save([path_output '/SPM.mat'], 'SPM');

         tmp = {};
         for f = 1:no_par
            tmp{f} = [path_output '/' sprintf('beta_%04d.nii', f)];
         end

         SPM.xY.P      = strvcat(tmp);
         SPM.xY.VY     = spm_vol(SPM.xY.P);
         SPM.xX.W      = eye(no_par);
         SPM.xX.K      = eye(no_par);
         SPM.xX.xKXs.X = eye(no_par);
         SPM.xX.iG     = [];
         SPM.xX.iB     = [];
         SPM.xX.Ic     = 0;

         save SPM SPM;

      end

   end

end

cd(path_manage);
