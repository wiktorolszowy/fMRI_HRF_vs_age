

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Figures for the analysis of the impact of age on HRFs and on the balloon parameters.
%%%%   Written by:  Wiktor Olszowy, University of Cambridge
%%%%   Contact:     wo222@cam.ac.uk
%%%%   Created:     October 2018 - December 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


paper             = 'HRF_age';
path_manage       = fgetl(fopen('path_manage.txt'));
path_scratch      = fgetl(fopen('path_scratch.txt'));
path_output       = [path_scratch '/analysis_output'];
dims_MNI          = [61 73 61];
slices            = [11 28 42 48];
subjects          = dir([path_scratch '/scans']);
subjects          = subjects(3:length(subjects));
subjects          = char({subjects.name});
subjects          = subjects(:, 1:8);
subjects          = unique(subjects, 'rows');
HRF_models        = cellstr(['canonical   '; 'canonical_TD'; 'FIR_32_05   '; 'FIR_16_1    '; 'FIR_24_1    ']);
HRF_models_t      = cellstr(['canonical           '; 'canonical+TD        '; 'FIR (32 x 0.5s bins)'; 'FIR (16 x 1s bins)  '; 'FIR (24 x 1s bins)  ']);
VOIs              = cellstr(['200'; '201'; '108'; '109'; '182'; '183'; '192'; '193']);
VOIs_t            = cellstr(['Right STG             ';
                             'Left STG              ';
                             'Right calcarine cortex';
                             'Left calcarine cortex ';
                             'Right precentral gyrus';
                             'Left precentral gyrus ';
                             'Right SMC             ';
                             'Left SMC              ']);
parameter_priors  = [0.65 0.41 0.98 0.32 0.34 -1 0];
parameter_names   = cellstr(['signal decay     '; 'autoregulation   '; 'transit time     '; 'Grubb''s exponent '; 'oxygen extraction'; 'intra:extra ratio'; 'neural efficacy  ']);
age_all_subjects  = textread('age_all_subjects.txt');
sensitivity_t     = cellstr(['transit_time     ';
                             'oxygen_extraction';
                             'all              ';
                             'precision        ']);
sensitivity_t2    = cellstr(['transit time 0.98->2    ';
                             'oxygen extr 0.34->0.4   ';
                             'all prior exp values x 2';
                             'prior precision / 10    ']);

cd(path_manage);
addpath(genpath([path_manage '/matlab_extra_functions']));
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');

load('combined_results/pos_rate.mat');
load('combined_results/avg_perc.mat');


%%%%%%%%%%%%%%%%%%%%%%%% 1st LEVEL: SPATIAL DISTRIBUTION OF SIGNIFICANT CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%
%-adjusted from the 'autocorr' study, as for other sizes/values, I encountered problems: some figures were a bit shifted...
for age   =  1:2
   if age == 1
      age_group = 'young';
   else
      age_group = 'old';
   end
   figure();
   for HRF_model_id = 1:length(HRF_models)
      %10>5 a trick, so that htitles are within the figure window and the subplots do not move
      subplot(10, 1, HRF_model_id+1);
      load(['combined_results/sp_dist_joint_F_' age_group '.mat']);
      sp_dist_joint_adjusted = zeros([4*dims_MNI(1) dims_MNI(2)]);
      for slice_id  = 1:length(slices)
         sp_dist_joint_adjusted((slice_id-1)*dims_MNI(1)+1:slice_id*dims_MNI(1), :) = sp_dist_joint_F((HRF_model_id-1)*dims_MNI(1)+1:HRF_model_id*dims_MNI(1), :, slices(slice_id));
      end
      imagesc(transpose(sp_dist_joint_adjusted), [0 60]);
      colormap gray; axis on; axis image; c = colorbar('FontSize', 6);
      set(gca, 'xticklabel', [], 'yticklabel', [], 'xtick', [], 'ytick', []);
      %-without latex, title not exactly in the middle...
      title(['\textbf{' HRF_models_t{HRF_model_id} '}'], 'FontSize', 6, 'interpreter', 'latex');
      if HRF_model_id == 1
         %-location can be estimated after getting position of 'title'
         fig_title = [age_group ' subjects'];
         if age == 1
            text(129,   -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 8.5);
         else
            text(127.2, -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 8.5);
         end
      end
      %-'Position': [left bottom width height]
      pos        = get(gca, 'Position');
      pos_bar    = get(c, 'Position');
      pos_bar(1) = 0.75;
      %-2nd term changed, because of pos(2) change, 3rd and 4th terms fixed as otherwise the size of bars can vary...
      pos_bar    = [pos_bar(1) pos(2)-0.023+0.1 0.0224 0.0595];
      set(c,   'Position', pos_bar);
      set(gca, 'Position', [0.04 pos(2) 0.66 0.21]);
   end
   set(gcf, 'units', 'inches', 'position', [0 0 4.0 ceil(10/7 * 7.5)]);
   figname = [paper '_1st_level_sp_dist_' age_group];
   print (['figures/' figname], '-dpdf');
   system(['pdfcrop figures/' figname '.pdf figures/' figname '.pdf']);
end


%%%%%%%%%%%%%%%%%%%%%%%% 2nd LEVEL: SPATIAL DISTRIBUTION OF SIGNIFICANT CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%
for age   =  1:2
   if age == 1
      age_group = 'young';
   else
      age_group = 'old';
   end
   figure();
   for HRF_model_id  = 1:length(HRF_models)
      HRF_model      = HRF_models{HRF_model_id};
      %10>5 a trick, so that htitles are within the figure window and the subplots do not move
      brain_template = rot90(niftiread('brain_parcellation/brain_template_MNI_3mm.nii'), 2);
      brain_template = brain_template / max(brain_template(:));
      cluster_binary_location = [path_output '/HRF_' HRF_model '/group_analysis_' age_group '/cluster_binary.nii'];
      if exist(cluster_binary_location, 'file') == 2
         cluster_binary = rot90(niftiread(cluster_binary_location), 2);
      else
         disp(['no cluster_binary at ' cluster_binary_location]);
         cluster_binary = zeros(dims_MNI);
      end
      sp_dist_joint_brain_template = zeros([4*dims_MNI(1) dims_MNI(2)]);
      sp_dist_joint_cluster_binary = zeros([4*dims_MNI(1) dims_MNI(2)]);
      for slice_id   = 1:length(slices)
         sp_dist_joint_brain_template((slice_id-1)*dims_MNI(1)+1:slice_id*dims_MNI(1), :) = brain_template(:, :, slices(slice_id));
         sp_dist_joint_cluster_binary((slice_id-1)*dims_MNI(1)+1:slice_id*dims_MNI(1), :) = cluster_binary(:, :, slices(slice_id));
      end
      sp_dist_joint_brain_template(sp_dist_joint_cluster_binary==1) = 1.2;
      subplot(10, 1, HRF_model_id+1);
      imagesc(transpose(sp_dist_joint_brain_template), [0 1.2]);
      color_palette        = gray();
      %-adding yellow for sig clusters
      color_palette(65, :) = [1 1 0];
      colormap(color_palette);
      axis on; axis image;
      set(gca, 'xticklabel', [], 'yticklabel', [], 'xtick', [], 'ytick', []);
      %-without latex, title not exactly in the middle...
      title(['\textbf{' HRF_models_t{HRF_model_id} '}'], 'FontSize', 6, 'interpreter', 'latex');
      if HRF_model_id == 1
         %-location can be estimated after getting position of 'title'
         fig_title = [age_group ' subjects'];
         if age == 1
            text(129,   -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 8.5);
         else
            text(127.2, -36, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 8.5);
         end
      end
      %-'Position': [left bottom width height]
      pos     = get(gca, 'Position');
      set(gca, 'Position', [0.04 pos(2) 0.66 0.21]);
      set(gcf, 'units', 'inches', 'position', [0 0 4.0 ceil(10/7 * 7.5)]);
   end
   figname = [paper '_2nd_level_sp_dist_' age_group];
   print (['figures/' figname], '-dpdf');
   system(['pdfcrop figures/' figname '.pdf figures/' figname '.pdf']);
end


%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATED PARAMETERS OF THE BALLOON MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
balloon_parameters = NaN([length(VOIs) length(parameter_priors) length(subjects)]);
HRF_model          = 'canonical_TD';
for i  = 1:length(VOIs)
   VOI = VOIs{i};
   load([path_scratch '/analysis_output/HRF_' HRF_model '/balloon_' VOI '.mat']);
   if length(Ep) ~= length(subjects)
      disp(['!!!!!!!!! ' VOI ' ' HRF_model ' ' num2str(length(Ep))]);
   end
   for parameter_id  = 1:length(parameter_priors)
      %-in almost all cases, length(Ep) = length(subjects)
      for subject_id = 1:length(Ep)
         %-checking if the balloon model was estimated, otherwise: NaN
         if ~isempty(Ep{subject_id})
            balloon_parameters(i, parameter_id, subject_id) = Ep{subject_id}(parameter_id);
         end
      end
   end
end
for part   =  1:2
   if part == 1
      VOIs_range = 1:4;
   else
      VOIs_range = 5:8;
   end
   figure('rend', 'painters', 'pos', [0 0 450 466], 'Visible', 'off');
   study_id_plot = 0;
   for parameter_id    = 1:length(parameter_priors)
      for i      = VOIs_range
         VOI_t   = VOIs_t{i};
         study_id_plot = study_id_plot + 1;
         subplot(7, 4, study_id_plot);
         axis off;
         %-[left bottom width height]
         gca_pos = get(gca, 'Position');
         %-lowering vertical spacing between subplots
         if parameter_id > 1
            gca_pos(2) = gca_pos(2) + (parameter_id-1)*0.004;
         end
         %-lowering horizontal spacing between subplots
         %-if sth changes in the figure, 0.043 might need to be adjusted...
         if i~=1 && i~=5
            gca_pos(1) = gca_pos(1) - rem(i-1, 4)*0.043;
         end
         axes('Position', gca_pos, 'Visible', 'off');
         par_for_subplot = balloon_parameters(i, parameter_id, :);
         par_for_subplot = par_for_subplot(:);
         idx     = ~isnan(par_for_subplot);
         %-when using 'plot', circles are open...
         scatter(age_all_subjects(idx), par_for_subplot(idx), 2.5, 'filled', 'o'); hold on;
         %-without 'box on', only x and y axes are shown (no closed box)
         box on;
         prior   = parameter_priors(parameter_id);
         par_all = balloon_parameters(:, parameter_id, :);
         axis([17 89 min([prior min(par_all(:))]) max([prior max(par_all(:))])]);
         lin_reg = fitlm(age_all_subjects(idx), par_for_subplot(idx));
         coefs   = lin_reg.Coefficients.Estimate;
         plot([0 100], [coefs(1) coefs(1)+100*coefs(2)], 'LineWidth', 1.25);
         [corr_R, corr_p] = corr(age_all_subjects(idx), par_for_subplot(idx));
         set(gca, 'FontSize', 4);
         if parameter_id < 7
            set(gca, 'XTickLabel', [], 'FontSize', 4);
         end
         if i~=1 && i~=5
            set(gca, 'YTickLabel', [], 'FontSize', 4);
         end
         title(['R=' sprintf('%.3f', corr_R) ', p=' sprintf('%.3f', corr_p)], 'FontSize', 5, 'FontWeight', 'normal');
         if i == 2 || i == 6
            text(1.024, 1.25, parameter_names{parameter_id}, 'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
         end
         if parameter_id == 1
            text(0.5,   1.4,  VOI_t,                         'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
         end
         %-SPM's prior expected value (horizontal line)
         plot([0 100], [prior prior], '--', 'Color', [0 0 0], 'LineWidth', 1.25); hold on;
         %-x-axis title: only for the last parameter -> last row of subplots
         if parameter_id == 7
            xlabel('Age');
         end
      end
   end
   figname = [paper '_balloon_parameters_' num2str(part)];
   print_to_svg_to_pdf(figname, path_manage);
end


%%%%%%%%%%%%%%%%%%%%%%%% BALLOON MODEL: SENSITIVITY ANALYSIS: HRF MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Left superior temporal gyrus
VOI           = VOIs{2};
balloon_parameters = NaN([6 length(parameter_priors) length(subjects)]);
figure('rend', 'painters', 'pos', [0 0 450 466], 'Visible', 'off');
study_id_plot = 0;
for HRF_model_id = 1:4
   HRF_model  = HRF_models{HRF_model_id};
   load([path_scratch '/analysis_output/HRF_' HRF_model '/balloon_' VOI '_sensitivity_HRF_models.mat']);
   if length(Ep) ~= length(subjects)
      disp(['!!!!!!!!! ' VOI ' ' HRF_model ' ' num2str(length(Ep))]);
   end
   for parameter_id  = 1:length(parameter_priors)
      %-in almost all cases, length(Ep) = length(subjects)
      for subject_id = 1:length(Ep)
         %-checking if the balloon model was estimated, otherwise: NaN
         if ~isempty(Ep{subject_id})
            balloon_parameters(HRF_model_id, parameter_id, subject_id) = Ep{subject_id}(parameter_id);
         end
      end
   end
end
for parameter_id    = 1:length(parameter_priors)
   for HRF_model_id = 1:4
      study_id_plot = study_id_plot + 1;
      subplot(7, 4, study_id_plot);
      axis off;
      %-[left bottom width height]
      gca_pos = get(gca, 'Position');
      %-lowering vertical spacing between subplots
      if parameter_id > 1
         gca_pos(2) = gca_pos(2) + (parameter_id-1)*0.004;
      end
      %-lowering horizontal spacing between subplots
      %-if sth changes in the figure, 0.043 might need to be adjusted...
      gca_pos(1)    = gca_pos(1) - (HRF_model_id-1)*0.043;
      axes('Position', gca_pos, 'Visible', 'off');
      par_for_subplot = balloon_parameters(HRF_model_id, parameter_id, :);
      par_for_subplot = par_for_subplot(:);
      idx     = ~isnan(par_for_subplot);
      %-when using 'plot', circles are open...
      scatter(age_all_subjects(idx), par_for_subplot(idx), 2.5, 'filled', 'o'); hold on;
      %-without 'box on', only x and y axes are shown (no closed box)
      box on;
      prior   = parameter_priors(parameter_id);
      par_all = balloon_parameters(:, parameter_id, :);
      axis([17 89 min([prior min(par_all(:))]) max([prior max(par_all(:))])]);
      lin_reg = fitlm(age_all_subjects(idx), par_for_subplot(idx));
      coefs   = lin_reg.Coefficients.Estimate;
      plot([0 100], [coefs(1) coefs(1)+100*coefs(2)], 'LineWidth', 1.25);
      [corr_R, corr_p] = corr(age_all_subjects(idx), par_for_subplot(idx));
      set(gca, 'FontSize', 4);
      if parameter_id < 7
         set(gca, 'XTickLabel', [], 'FontSize', 4);
      end
      if HRF_model_id > 1
         set(gca, 'YTickLabel', [], 'FontSize', 4);
      end
      title(['R=' sprintf('%.3f', corr_R) ', p=' sprintf('%.3f', corr_p)], 'FontSize', 5, 'FontWeight', 'normal');
      if HRF_model_id == 2
         text(1.024, 1.25, parameter_names{parameter_id}, 'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
      end
      if parameter_id == 1
         text(0.5,   1.4,  HRF_models_t{HRF_model_id},    'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
      end
      %-SPM's prior expected value (horizontal line)
      plot([0 100], [prior prior], '--', 'Color', [0 0 0], 'LineWidth', 1.25); hold on;
      %-x-axis title: only for the last parameter -> last row of subplots
      if parameter_id == 7
         xlabel('Age');
      end
   end
end
figname = [paper '_balloon_parameters_sensitivity_HRF_models'];
print_to_svg_to_pdf(figname, path_manage);


%%%%%%%%%%%%%%%%%%%%%%%% BALLOON MODEL: SENSITIVITY ANALYSIS: PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HRF_model     = 'canonical_TD';
%-Left superior temporal gyrus
VOI           = VOIs{2};
balloon_parameters = NaN([length(sensitivity_t) length(parameter_priors) length(subjects)]);
figure('rend', 'painters', 'pos', [0 0 450 466], 'Visible', 'off');
study_id_plot = 0;
for sensitivity_case = 1:length(sensitivity_t)
   load([path_scratch '/analysis_output/HRF_' HRF_model '/balloon_' VOI '_sensitivity_' sensitivity_t{sensitivity_case} '.mat']);
   if length(Ep) ~= length(subjects)
      disp(['!!!!!!!!! ' VOI ' ' HRF_model ' ' num2str(length(Ep))]);
   end
   for parameter_id  = 1:length(parameter_priors)
      %-in almost all cases, length(Ep) = length(subjects)
      for subject_id = 1:length(Ep)
         %-checking if the balloon model was estimated, otherwise: NaN
         if ~isempty(Ep{subject_id})
            balloon_parameters(sensitivity_case, parameter_id, subject_id) = Ep{subject_id}(parameter_id);
         end
      end
   end
end
for parameter_id = 1:length(parameter_priors)
   for sensitivity_case          =  1:length(sensitivity_t)
      parameter_priors_aux       =  parameter_priors;
      if sensitivity_case        == 1
         parameter_priors_aux(3) =  2;
      elseif sensitivity_case    == 2
         parameter_priors_aux(5) =  0.4;
      elseif sensitivity_case    == 3
         parameter_priors_aux    =  2*parameter_priors_aux;
      end
      study_id_plot = study_id_plot + 1;
      subplot(7, 4, study_id_plot);
      axis off;
      %-[left bottom width height]
      gca_pos = get(gca, 'Position');
      %-lowering vertical spacing between subplots
      if parameter_id > 1
         gca_pos(2) = gca_pos(2) + (parameter_id-1)*0.004;
      end
      %-lowering horizontal spacing between subplots
      %-if sth changes in the figure, 0.043 might need to be adjusted...
      gca_pos(1)      = gca_pos(1) - (sensitivity_case-1)*0.043;
      axes('Position', gca_pos, 'Visible', 'off');
      par_for_subplot = balloon_parameters(sensitivity_case, parameter_id, :);
      par_for_subplot = par_for_subplot(:);
      idx     = ~isnan(par_for_subplot);
      %-when using 'plot', circles are open...
      scatter(age_all_subjects(idx), par_for_subplot(idx), 2.5, 'filled', 'o'); hold on;
      %-without 'box on', only x and y axes are shown (no closed box)
      box on;
      prior   = parameter_priors_aux(parameter_id);
      par_all = balloon_parameters(:, parameter_id, :);
      axis([17 89 min([prior min(par_all(:))]) max([prior max(par_all(:))])]);
      lin_reg = fitlm(age_all_subjects(idx), par_for_subplot(idx));
      coefs   = lin_reg.Coefficients.Estimate;
      plot([0 100], [coefs(1) coefs(1)+100*coefs(2)], 'LineWidth', 1.25);
      [corr_R, corr_p] = corr(age_all_subjects(idx), par_for_subplot(idx));
      set(gca, 'FontSize', 4);
      if parameter_id < 7
         set(gca, 'XTickLabel', [], 'FontSize', 4);
      end
      if sensitivity_case > 1
         set(gca, 'YTickLabel', [], 'FontSize', 4);
      end
      title(['R=' sprintf('%.3f', corr_R) ', p=' sprintf('%.3f', corr_p)], 'FontSize', 5, 'FontWeight', 'normal');
      if sensitivity_case == 2
         text(1.024, 1.25, parameter_names{parameter_id},    'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
      end
      if parameter_id == 1
         text(0.5,   1.4,  sensitivity_t2{sensitivity_case}, 'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
      end
      %-SPM's prior expected value (horizontal line)
      plot([0 100], [prior prior], '--', 'Color', [0 0 0], 'LineWidth', 1.25); hold on;
      %-x-axis title: only for the last parameter -> last row of subplots
      if parameter_id == 7
         xlabel('Age');
      end
   end
end
figname = [paper '_balloon_parameters_sensitivity_priors'];
print_to_svg_to_pdf(figname, path_manage);

cd(path_manage);
