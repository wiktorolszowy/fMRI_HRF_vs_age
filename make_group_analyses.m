

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Group analyses for CamCAN data.
%%%%   Written by:    Wiktor Olszowy, University of Cambridge
%%%%   Contact:       wo222@cam.ac.uk
%%%%   Created:       October 2018 - December 2018
%%%%   Adapted from:  https://github.com/wanderine/ParametricMultisubjectfMRI/blob/master/SPM/run_random_group_analyses_onesamplettest_1.m
%%%%                  http://www.fil.ion.ucl.ac.uk/spm/data/face_rfx/
%%%%                  https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/make_group_analyses_random_effects.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage    = fgetl(fopen('path_manage.txt'));
path_scratch   = fgetl(fopen('path_scratch.txt'));
HRF_models     = cellstr(['canonical   '; 'canonical_TD'; 'FIR_32_05   '; 'FIR_16_1    '; 'FIR_24_1    ']);

addpath(genpath([path_manage '/matlab_extra_functions']));
addpath('/applications/spm/spm12_7219');

spm('Defaults', 'fMRI');
spm_jobman('initcfg');


%-run in R to choose 30 young and 30 old subjects:
%{
path_scratch   = readLines("path_scratch.txt")
subjects       = list.files(paste0(path_scratch, "/scans"))
subjects       = unique(substr(subjects, 1, 8))
A              = read.csv("CamCAN_subjects_info.csv")
A              = A[which(A$CCID %in% subjects), ]
subjects_young = A[order(A$Age),]$CCID[1:30]
subjects_old   = A[order(A$Age),]$CCID[(length(subjects)-30+1):length(subjects)]
subjects_young = paste0("'", subjects_young, "'; ", collapse='')
subjects_old   = paste0("'", subjects_old,   "'; ", collapse='')
age            = A$Age[which(A$CCID%in%subjects)]
write(age, "age_all_subjects.txt", ncolumns=1)
%}
subjects_young = cellstr(['CC110037'; 'CC110182'; 'CC120376'; 'CC120409'; 'CC120462'; 'CC121111'; 'CC120061'; 'CC120123'; 'CC120550'; 'CC110606'; 'CC120987'; 'CC121685'; 'CC120286'; 'CC120347'; 'CC110056'; 'CC110126'; 'CC110098'; 'CC110101'; 'CC120276'; 'CC120727'; 'CC120816'; 'CC110033'; 'CC120208'; 'CC120234'; 'CC120795'; 'CC121194'; 'CC122620'; 'CC110174'; 'CC110187'; 'CC110411']);
subjects_old   = cellstr(['CC721704'; 'CC710088'; 'CC710154'; 'CC710342'; 'CC710382'; 'CC710548'; 'CC710566'; 'CC720358'; 'CC720986'; 'CC721891'; 'CC710679'; 'CC720290'; 'CC720516'; 'CC721434'; 'CC722891'; 'CC710099'; 'CC710131'; 'CC710446'; 'CC710551'; 'CC710591'; 'CC711244'; 'CC711245'; 'CC720400'; 'CC721374'; 'CC722216'; 'CC723395'; 'CC712027'; 'CC721224'; 'CC721532'; 'CC711035']);

HRF_model      = HRF_models{HRF_model_id};
path_output    = [path_scratch '/analysis_output/HRF_' HRF_model];

cd(path_output);

if age == 1
   age_group   = 'young';
else
   age_group   = 'old';
end

SPM_mat_location = cellstr(fullfile(path_output, ['group_analysis_' age_group], 'SPM.mat'));

%-when the folder is not empty, problems occur
%-when repeating a run, slight differences in results, probably some randomization involved (?)

system('rm -rf group_analysis');
system(['rm -rf group_analysis_' age_group]);
system(['mkdir group_analysis_' age_group]);

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

clear jobs;

jobs{1}.stats{1}.factorial_design.dir                              = cellstr(fullfile(path_output, ['group_analysis_' age_group]));
%-it seems there are errors when the name is not specified
jobs{1}.stats{1}.factorial_design.des.fd.fact.name                 = 'Basis';
jobs{1}.stats{1}.factorial_design.des.fd.fact.levels               = no_par;
jobs{1}.stats{1}.factorial_design.des.fd.fact.dept                 = 1;
for par_id = 1:no_par
   jobs{1}.stats{1}.factorial_design.des.fd.icell(par_id).levels   = par_id;
   if par_id < 10
      coef_maps                                                    = cellstr(spm_select('FPListRec', fullfile(path_output), ['^beta_000' num2str(par_id) '.nii']));
   else
      coef_maps                                                    = cellstr(spm_select('FPListRec', fullfile(path_output), ['^beta_00'  num2str(par_id) '.nii']));
   end
   if strcmp(age_group, 'young')
      jobs{1}.stats{1}.factorial_design.des.fd.icell(par_id).scans = coef_maps(find(contains(coef_maps, subjects_young)));
   else
      jobs{1}.stats{1}.factorial_design.des.fd.icell(par_id).scans = coef_maps(find(contains(coef_maps, subjects_old)));
   end
end
jobs{1}.stats{2}.fmri_est.spmmat                                   = SPM_mat_location;
jobs{1}.stats{3}.con.spmmat                                        = SPM_mat_location;
jobs{1}.stats{3}.con.consess{1}.fcon.name                          = 'F-contrast';
jobs{1}.stats{3}.con.consess{1}.fcon.weights                       = eye(no_par);    %for two samples: -[1 -1 0 0];
jobs{1}.stats{4}.results.spmmat                                    = SPM_mat_location;
jobs{1}.stats{4}.results.conspec.contrasts                         = Inf;
jobs{1}.stats{4}.results.conspec.threshdesc                        = 'FWE';

spm_jobman('run', jobs);

%-get F-map
V           = spm_vol([path_output '/group_analysis_' age_group '/spmF_0001.nii']);
[tmap,aa]   = spm_read_vols(V);

%-calculate cluster extent threshold
[k,Pc]      = corrclusth(SPM, 0.001, 0.05, 1:100000);
df          = [no_par SPM.xX.erdf];
u           = spm_u(0.001, df, 'F');
indices     = find(tmap>u);

%-get the size of the largest cluster
max_cluster = max_extent(tmap, indices);

if max_cluster >= k
   indices     = find(tmap>u);
   disp([age_group ' ' HRF_model ' : significant!']);
else
   indices     = '';
end
fid            = fopen('indices.txt', 'wt');
fprintf(fid, '%8.0f\n', indices);

mask                    = niftiread('mask.nii');
cluster_binary          = zeros(size(mask));
indices                 = fscanf(fopen('indices.txt', 'r'), '%d');
cluster_binary(indices) = 1;
niftiwrite(cluster_binary, 'cluster_binary.nii');

cd(path_manage)
