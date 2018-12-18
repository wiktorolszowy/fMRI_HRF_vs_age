

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Estimating the balloon model on already prepared VOIs.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       November 2018 - December 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage  = fgetl(fopen('path_manage.txt'));
path_scratch = fgetl(fopen('path_scratch.txt'));
HRF_models   = cellstr(['canonical   '; 'canonical_TD'; 'FIR_32_05   '; 'FIR_16_1    '; 'FIR_24_1    ']);
subjects     = dir([path_scratch '/scans']);
subjects     = subjects(3:length(subjects));
subjects     = char({subjects.name});
subjects     = subjects(:, 1:8);
subjects     = unique(subjects, 'rows');
p            = parpool(12);
VOIs         = cellstr(['200'; '201'; '108'; '109'; '182'; '183'; '192'; '193']);
HRF_model    = 'canonical_TD';

%-VOI numbers refer to /applications/spm/spm12_7219/tpm/labels_Neuromorphometrics.nii

addpath('/applications/spm/spm12_7219');
addpath(genpath([path_manage '/matlab_extra_functions']));

%-initialise SPM defaults
spm('Defaults', 'fMRI');
spm_jobman('initcfg');

path_output  = [path_scratch '/analysis_output/HRF_' HRF_model];

cd(path_output);


for i  = 1:length(VOIs)

   VOI = VOIs{i};
   disp(VOI);

   Ep  = {};
   Cp  = {};
   K1  = {};
   K2  = {};
   M0  = {};
   M1  = {};

   parfor subject_id = 1:length(subjects)

      subject  = subjects(subject_id, :);

      disp(subject_id);

      SPM_file = [path_output '/' subject '/VOI_' VOI '/SPM.mat'];
      VOI_file = [path_output '/' subject '/VOI_' VOI '/VOI_' VOI '_1.mat'];

      %-if the VOI was not produced, the 'subject_id' location in the output vector will be empty
      %-importantly, if this is the last 'subject_id', the output vector will be shorter!!!
      if exist(SPM_file, 'file') && exist(VOI_file, 'file')

         [a, b, c, d, e, f] = rik_hdm(SPM_file, VOI_file);

         Ep{subject_id} = a;
         Cp{subject_id} = b;
         K1{subject_id} = c;
         K2{subject_id} = d;
         M0{subject_id} = e;
         M1{subject_id} = f;

      end

   end

   save(['balloon_' VOI], 'Ep', 'Cp', 'K1', 'K2', 'M0', 'M1');

end

cd(path_manage);
