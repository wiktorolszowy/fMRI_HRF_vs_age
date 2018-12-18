

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Combining results from single runs to 'mat' files. Necessary for making figures.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       October 2018 - December 2018
%%%%   Adapted from:  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/make_combined_results_from_single_runs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage       = fgetl(fopen('path_manage.txt'));
path_scratch      = fgetl(fopen('path_scratch.txt'));
path_output       = [path_scratch '/analysis_output'];
subjects          = dir([path_scratch '/scans']);
subjects          = subjects(3:length(subjects));
subjects          = char({subjects.name});
subjects          = subjects(:, 1:8);
subjects          = unique(subjects, 'rows');
no_subjects       = length(subjects);
HRF_models        = cellstr(['canonical   '; 'canonical_TD'; 'FIR_32_05   '; 'FIR_16_1    '; 'FIR_24_1    ']);
combined          = (-1)*ones(no_subjects, length(HRF_models));
combined_perc     = (-1)*ones(no_subjects, length(HRF_models));
dims              = size(combined);
pos_rate          = (-1)*ones(dims(2), 1);
avg_perc          = (-1)*ones(dims(2), 1);
range_HRF_models  = 1:length(HRF_models);
dims_MNI          = [61 73 61];
p                 = parpool(12);


parfor subject_id = 1:no_subjects

   subject        = subjects(subject_id, :);

   combined_parfor      = (-1)*ones(length(HRF_models), 1);
   combined_perc_parfor = (-1)*ones(length(HRF_models), 1);

   for HRF_model_id = range_HRF_models

      HRF_model     = HRF_models{HRF_model_id};

      cd([path_output '/HRF_' HRF_model]);

      if exist([subject '/cluster_binary_MNI.nii'], 'file') == 2
         cd(subject);
      else
         disp([pwd ' ' subject]);
         continue
      end

      MNI_sig                       = niftiread('cluster_binary_MNI.nii');
      MNI_mask                      = niftiread('mask.nii');
      combined_parfor(HRF_model_id) = sum(MNI_sig(:));

      if sum(MNI_mask(:)) > 0

         combined_perc_parfor(HRF_model_id) = sum(MNI_sig(:)) / sum(MNI_mask(:));

      else

         disp(pwd);

      end

   end

   combined     (subject_id, :) = combined_parfor;
   combined_perc(subject_id, :) = combined_perc_parfor;

end

for i = 1:dims(2)

   over_sub       = combined     (:, i);
   over_sub_perc  = combined_perc(:, i);

   %-(-0.5) chosen for numerical reasons; in fact, we only want to distinguish >=0 from <0
   if sum(over_sub>-0.5) > 0
      pos_rate(i) =  sum(over_sub>0)/sum(over_sub>-0.5);
      %-sum(over_sub<-0.5) considered, as that many times the default (-1) is subtracted; (-1) appears in 'combined' for non-subjects
      avg_perc(i) = (sum(over_sub_perc) + sum(over_sub<-0.5)) / sum(over_sub>-0.5);
   end

end

cd(path_manage);
save('combined_results/combined',      'combined');
save('combined_results/combined_perc', 'combined_perc');
save('combined_results/pos_rate',      'pos_rate');
save('combined_results/avg_perc',      'avg_perc');


%-saving 'mat' files which show spatial distribution of significant clusters (from 1st level analyses)

subjects_young  = cellstr(['CC110037'; 'CC110182'; 'CC120376'; 'CC120409'; 'CC120462'; 'CC121111'; 'CC120061'; 'CC120123'; 'CC120550'; 'CC110606'; 'CC120987'; 'CC121685'; 'CC120286'; 'CC120347'; 'CC110056'; 'CC110126'; 'CC110098'; 'CC110101'; 'CC120276'; 'CC120727'; 'CC120816'; 'CC110033'; 'CC120208'; 'CC120234'; 'CC120795'; 'CC121194'; 'CC122620'; 'CC110174'; 'CC110187'; 'CC110411']);
subjects_old    = cellstr(['CC721704'; 'CC710088'; 'CC710154'; 'CC710342'; 'CC710382'; 'CC710548'; 'CC710566'; 'CC720358'; 'CC720986'; 'CC721891'; 'CC710679'; 'CC720290'; 'CC720516'; 'CC721434'; 'CC722891'; 'CC710099'; 'CC710131'; 'CC710446'; 'CC710551'; 'CC710591'; 'CC711244'; 'CC711245'; 'CC720400'; 'CC721374'; 'CC722216'; 'CC723395'; 'CC712027'; 'CC721224'; 'CC721532'; 'CC711035']);
sp_dist_joint_F = zeros(length(HRF_models)*dims_MNI(1), dims_MNI(2), dims_MNI(3));

for age   =  1:2

   if age == 1
      age_group = 'young';
   else
      age_group = 'old';
   end

   for HRF_model_id     = 1:length(HRF_models)

      HRF_model         = HRF_models{HRF_model_id};
      sp_dist_F         = zeros(dims_MNI);
      no_without_output = 0;

      for subject_id    = 1:30

         if strcmp(age_group, 'young')
            subject                    = subjects_young{subject_id};
         else
            subject                    = subjects_old{subject_id};
         end

         cluster_binary_F_MNI_location = [path_output '/HRF_' HRF_model '/' subject '/cluster_binary_MNI.nii'];

         if exist(cluster_binary_F_MNI_location, 'file') == 2
            cluster_binary_F_MNI       = niftiread(cluster_binary_F_MNI_location);
            sp_dist_F                  = sp_dist_F + cluster_binary_F_MNI;
         else
            disp(['problems with sp_dists_F for ' cluster_binary_F_MNI_location]);
            no_without_output          = no_without_output + 1;
         end

      end

      sp_dist_F = rot90(sp_dist_F, 2);
      sp_dist_joint_F((HRF_model_id-1)*dims_MNI(1)+1:HRF_model_id*dims_MNI(1), :, :) = 100 * sp_dist_F / (30 - no_without_output);

   end

   save(['combined_results/sp_dist_joint_F_' age_group '.mat'], 'sp_dist_joint_F');

end
