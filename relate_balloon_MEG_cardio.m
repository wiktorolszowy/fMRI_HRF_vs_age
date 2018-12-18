

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Relating balloon parameters to MEG delays and cardio data.
%%%%   Written by:  Wiktor Olszowy, University of Cambridge
%%%%   Contact:     wo222@cam.ac.uk
%%%%   Created:     November 2018 - December 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage      = fgetl(fopen('path_manage.txt'));
path_scratch     = fgetl(fopen('path_scratch.txt'));
VOIs             = cellstr(['200'; '201'; '108'; '109'; '182'; '183'; '192'; '193']);
pars_balloon_t   = cellstr(['signal decay     '; 'autoregulation   '; 'transit time     '; 'Grubb''s exponent '; 'oxygen extraction'; 'intra:extra ratio'; 'neural efficacy  ']);
pars_MEG_t       = cellstr(['MEG constant delay [ms]  '; 'MEG cumulative delay [ms]'; 'MEG amplitude            ']);
pars_cardio_t    = cellstr(['Mean systolic [mmHg]   '; 'Mean diastolic [mmHg]  '; 'Mean pulse [mmHg]      '; 'Pulse pressure [mmHg]  ']);
paper            = 'HRF_age';

addpath(genpath([path_manage '/matlab_extra_functions']));


%-defining subjects for which balloon model was estimated, and loading these estimates
subjects_balloon = dir([path_scratch '/scans']);
subjects_balloon = subjects_balloon(3:length(subjects_balloon));
subjects_balloon = char({subjects_balloon.name});
subjects_balloon = subjects_balloon(:, 1:8);
subjects_balloon = cellstr(unique(subjects_balloon, 'rows'));
combined_balloon = NaN([length(VOIs) length(pars_balloon_t) length(subjects_balloon)]);
HRF_model        = 'canonical_TD';
for i  = 1:length(VOIs)
   VOI = VOIs{i};
   load([path_scratch '/analysis_output/HRF_' HRF_model '/balloon_' VOI '.mat']);
   if length(Ep) ~= length(subjects_balloon)
      disp(['!!!!!!!!! ' VOI ' ' HRF_model ' ' num2str(length(Ep))]);
   end
   for par_balloon_id = 1:length(pars_balloon_t)
      %-in almost all cases, length(Ep) = length(subjects_balloon)
      for subject_id  = 1:length(Ep)
         %-checking if the balloon model was estimated, otherwise: NaN
         if ~isempty(Ep{subject_id})
            combined_balloon(i, par_balloon_id, subject_id) = Ep{subject_id}(par_balloon_id);
         end
      end
   end
end


%-defining subjects for which MEG delays were estimated, and loading these estimates
combined_MEG     = load('MEG_cardio/avpass_latency_amp_table');
names_MEG        = combined_MEG.d.Properties.VariableNames;
combined_MEG     = table2array(combined_MEG.d)';
subjects_MEG     = readtable('CamCAN_subjects_info.csv');
subjects_MEG     = subjects_MEG.CCID;


%-defining subjects for which cardio data were obtained, and loading these data
load('MEG_cardio/cc700_cardio.mat');
combined_cardio  = cardio;
names_cardio     = cardio.Properties.VariableNames;
subjects_cardio  = combined_cardio.CCID;
combined_cardio  = [combined_cardio.bp_sys_mean'; combined_cardio.bp_dia_mean'; combined_cardio.pulse_mean'; combined_cardio.bp_pp'];


%-identifying balloon outliers, following Price ea. 2017 Nature Comm:
%-2.5 chosen rather than 1.5 as more dimensions are considered (7 rather than 2), the value (2.5) chosen arbitrarily
%-for left STG
outliers_balloon_aud = [find(combined_balloon(2, 1, :) < prctile(combined_balloon(2, 1, :), 25) - 2.5*iqr(combined_balloon(2, 1, :)));
                        find(combined_balloon(2, 1, :) > prctile(combined_balloon(2, 1, :), 25) + 2.5*iqr(combined_balloon(2, 1, :)));
                        find(combined_balloon(2, 2, :) < prctile(combined_balloon(2, 2, :), 25) - 2.5*iqr(combined_balloon(2, 2, :)));
                        find(combined_balloon(2, 2, :) > prctile(combined_balloon(2, 2, :), 25) + 2.5*iqr(combined_balloon(2, 2, :)));
                        find(combined_balloon(2, 3, :) < prctile(combined_balloon(2, 3, :), 25) - 2.5*iqr(combined_balloon(2, 3, :)));
                        find(combined_balloon(2, 3, :) > prctile(combined_balloon(2, 3, :), 25) + 2.5*iqr(combined_balloon(2, 3, :)));
                        find(combined_balloon(2, 4, :) < prctile(combined_balloon(2, 4, :), 25) - 2.5*iqr(combined_balloon(2, 4, :)));
                        find(combined_balloon(2, 4, :) > prctile(combined_balloon(2, 4, :), 25) + 2.5*iqr(combined_balloon(2, 4, :)));
                        find(combined_balloon(2, 5, :) < prctile(combined_balloon(2, 5, :), 25) - 2.5*iqr(combined_balloon(2, 5, :)));
                        find(combined_balloon(2, 5, :) > prctile(combined_balloon(2, 5, :), 25) + 2.5*iqr(combined_balloon(2, 5, :)));
                        find(combined_balloon(2, 6, :) < prctile(combined_balloon(2, 6, :), 25) - 2.5*iqr(combined_balloon(2, 6, :)));
                        find(combined_balloon(2, 6, :) > prctile(combined_balloon(2, 6, :), 25) + 2.5*iqr(combined_balloon(2, 6, :)));
                        find(combined_balloon(2, 7, :) < prctile(combined_balloon(2, 7, :), 25) - 2.5*iqr(combined_balloon(2, 7, :)));
                        find(combined_balloon(2, 7, :) > prctile(combined_balloon(2, 7, :), 25) + 2.5*iqr(combined_balloon(2, 7, :)));];
outliers_balloon_aud = unique(outliers_balloon_aud);
%-for left calcarine cortex
outliers_balloon_vis = [find(combined_balloon(4, 1, :) < prctile(combined_balloon(4, 1, :), 25) - 2.5*iqr(combined_balloon(4, 1, :)));
                        find(combined_balloon(4, 1, :) > prctile(combined_balloon(4, 1, :), 25) + 2.5*iqr(combined_balloon(4, 1, :)));
                        find(combined_balloon(4, 2, :) < prctile(combined_balloon(4, 2, :), 25) - 2.5*iqr(combined_balloon(4, 2, :)));
                        find(combined_balloon(4, 2, :) > prctile(combined_balloon(4, 2, :), 25) + 2.5*iqr(combined_balloon(4, 2, :)));
                        find(combined_balloon(4, 3, :) < prctile(combined_balloon(4, 3, :), 25) - 2.5*iqr(combined_balloon(4, 3, :)));
                        find(combined_balloon(4, 3, :) > prctile(combined_balloon(4, 3, :), 25) + 2.5*iqr(combined_balloon(4, 3, :)));
                        find(combined_balloon(4, 4, :) < prctile(combined_balloon(4, 4, :), 25) - 2.5*iqr(combined_balloon(4, 4, :)));
                        find(combined_balloon(4, 4, :) > prctile(combined_balloon(4, 4, :), 25) + 2.5*iqr(combined_balloon(4, 4, :)));
                        find(combined_balloon(4, 5, :) < prctile(combined_balloon(4, 5, :), 25) - 2.5*iqr(combined_balloon(4, 5, :)));
                        find(combined_balloon(4, 5, :) > prctile(combined_balloon(4, 5, :), 25) + 2.5*iqr(combined_balloon(4, 5, :)));
                        find(combined_balloon(4, 6, :) < prctile(combined_balloon(4, 6, :), 25) - 2.5*iqr(combined_balloon(4, 6, :)));
                        find(combined_balloon(4, 6, :) > prctile(combined_balloon(4, 6, :), 25) + 2.5*iqr(combined_balloon(4, 6, :)));
                        find(combined_balloon(4, 7, :) < prctile(combined_balloon(4, 7, :), 25) - 2.5*iqr(combined_balloon(4, 7, :)));
                        find(combined_balloon(4, 7, :) > prctile(combined_balloon(4, 7, :), 25) + 2.5*iqr(combined_balloon(4, 7, :)));];
outliers_balloon_vis = unique(outliers_balloon_vis);


%-identifying MEG outliers, following Price ea. 2017 Nature Comm:
%-outliers either in the constant or in the cumulative delay estimate
%-for auditory cortex
outliers_MEG_aud = [find(combined_MEG(1, :) < prctile(combined_MEG(1, :), 25) - 1.5*iqr(combined_MEG(1, :)))';
                    find(combined_MEG(1, :) > prctile(combined_MEG(1, :), 75) + 1.5*iqr(combined_MEG(1, :)))';
                    find(combined_MEG(2, :) < prctile(combined_MEG(2, :), 25) - 1.5*iqr(combined_MEG(2, :)))';
                    find(combined_MEG(2, :) > prctile(combined_MEG(2, :), 75) + 1.5*iqr(combined_MEG(2, :)))'];
outliers_MEG_aud = unique(outliers_MEG_aud);
%-for visual cortex
outliers_MEG_vis = [find(combined_MEG(4, :) < prctile(combined_MEG(4, :), 25) - 1.5*iqr(combined_MEG(4, :)))';
                    find(combined_MEG(4, :) > prctile(combined_MEG(4, :), 75) + 1.5*iqr(combined_MEG(4, :)))';
                    find(combined_MEG(5, :) < prctile(combined_MEG(5, :), 25) - 1.5*iqr(combined_MEG(5, :)))';
                    find(combined_MEG(5, :) > prctile(combined_MEG(5, :), 75) + 1.5*iqr(combined_MEG(5, :)))'];
outliers_MEG_vis = unique(outliers_MEG_vis);


%%%%%%%%%%%%%%%%%%%%%%%% BALLOON MODEL + MEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combined_balloon_MEG_aud = [];
combined_balloon_MEG_vis = [];

for subject_balloon_id = 1:length(combined_balloon)
   subject_MEG_id      = find(strcmp(subjects_balloon{subject_balloon_id}, subjects_MEG));
   %-in 'combined_balloon' NaNs always appear together (column-wise), in 'combined_MEG' too
   if ~isempty(subject_MEG_id) && ~isnan(combined_balloon(2, 1, subject_balloon_id)) && ~isnan(combined_MEG(1, subject_MEG_id))
      %-left STG
      if ~isnan(combined_balloon(2, 1, subject_balloon_id)) && ~ismember(subject_balloon_id, outliers_balloon_aud) && ~ismember(subject_MEG_id, outliers_MEG_aud)
         combined_balloon_MEG_aud(:, end+1) = [combined_balloon(2, 1:7, subject_balloon_id) combined_MEG(1:3, subject_MEG_id)'];
      end
   end
   if ~isempty(subject_MEG_id) && ~isnan(combined_balloon(4, 1, subject_balloon_id)) && ~isnan(combined_MEG(1, subject_MEG_id))
      %-left calcarine cortex
      if ~isnan(combined_balloon(4, 1, subject_balloon_id)) && ~ismember(subject_balloon_id, outliers_balloon_vis) && ~ismember(subject_MEG_id, outliers_MEG_vis)
         combined_balloon_MEG_vis(:, end+1) = [combined_balloon(4, 1:7, subject_balloon_id) combined_MEG(4:6, subject_MEG_id)'];
      end
   end
end

%-auditory cortex: correlations between balloon and MEG
disp(['   '; 'aud'; '   ']);
for par_balloon_id = 1:length(pars_balloon_t)
   for par_MEG_id  = 1:length(pars_MEG_t)
      [corr_R, corr_p] = corr(combined_balloon_MEG_aud(par_balloon_id, :)', combined_balloon_MEG_aud(7+par_MEG_id, :)');
      disp([pars_balloon_t{par_balloon_id} ' ' pars_MEG_t{par_MEG_id} ' ' num2str(corr_p)]);
   end
end

disp(['   '; 'vis'; '   ']);
%-visual cortex: correlations between balloon and MEG
for par_balloon_id = 1:length(pars_balloon_t)
   for par_MEG_id  = 1:length(pars_MEG_t)
      [corr_R, corr_p] = corr(combined_balloon_MEG_vis(par_balloon_id, :)', combined_balloon_MEG_vis(7+par_MEG_id, :)');
      disp([pars_balloon_t{par_balloon_id} ' ' pars_MEG_t{par_MEG_id} ' ' num2str(corr_p)]);
   end
end

%-making figures
for cortex_id   =  1:2
   if cortex_id == 1
      combined_balloon_MEG = combined_balloon_MEG_aud;
      cortex = 'aud';
   else
      combined_balloon_MEG = combined_balloon_MEG_vis;
      cortex = 'vis';
   end
   figure('rend', 'painters', 'pos', [0 0 450 466], 'Visible', 'off');
   study_id_plot = 0;
   for par_balloon_id = 1:length(pars_balloon_t)
      for par_MEG_id  = 1:length(pars_MEG_t)
         study_id_plot = study_id_plot + 1;
         subplot(7, 3, study_id_plot);
         x_val = combined_balloon_MEG(7+par_MEG_id,   :)';
         y_val = combined_balloon_MEG(par_balloon_id, :)';
         axis off;
         %-[left bottom width height]
         gca_pos = get(gca, 'Position');
         %-lowering vertical spacing between subplots
         if par_balloon_id > 1
            gca_pos(2) = gca_pos(2) + (par_balloon_id-1)*0.004;
         end
         %-lowering horizontal spacing between subplots
         %-if sth changes in the figure, 0.043 might need to be adjusted...
         gca_pos(1)      = gca_pos(1) - (par_MEG_id-1)*0.043;
         axes('Position', gca_pos, 'Visible', 'off');
         par_for_subplot = combined_balloon_MEG(par_balloon_id, par_MEG_id, :);
         par_for_subplot = par_for_subplot(:);
         %-when using 'plot', circles are open...
         scatter(x_val, y_val, 2.5, 'filled', 'o'); hold on;
         %-without 'box on', only x and y axes are shown (no closed box)
         box on;
         axis([min(combined_balloon_MEG(7+par_MEG_id, :)) max(combined_balloon_MEG(7+par_MEG_id, :)) min(combined_balloon_MEG(par_balloon_id, :)) max(combined_balloon_MEG(par_balloon_id, :))]);
         lin_reg = fitlm(x_val, y_val);
         coefs   = lin_reg.Coefficients.Estimate;
         plot([min(x_val) max(x_val)], [coefs(1)+min(x_val)*coefs(2) coefs(1)+max(x_val)*coefs(2)], 'LineWidth', 1.25);
         [corr_R, corr_p] = corr(x_val, y_val);
         set(gca, 'FontSize', 4);
         if par_balloon_id < 7
            set(gca, 'XTickLabel', [], 'FontSize', 4);
         end
         if par_MEG_id > 1
            set(gca, 'YTickLabel', [], 'FontSize', 4);
         end
         title(['R=' sprintf('%.3f', corr_R) ', p=' sprintf('%.3f', corr_p)], 'FontSize', 5, 'FontWeight', 'normal');
         if par_MEG_id == 2
            text(0.5, 1.25, pars_balloon_t{par_balloon_id}, 'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
         end
         %-x-axis title: only for the last parameter -> last row of subplots
         if par_balloon_id == 7
            xlabel(pars_MEG_t{par_MEG_id}, 'FontSize', 5, 'FontWeight', 'bold');
         end
      end
   end
   figname = [paper '_balloon_MEG_' cortex];
   print_to_svg_to_pdf(figname, path_manage);
end


%%%%%%%%%%%%%%%%%%%%%%%% BALLOON MODEL + cardio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combined_balloon_cardio_aud = [];
combined_balloon_cardio_vis = [];

for subject_balloon_id = 1:length(combined_balloon)
   subject_cardio_id   = find(strcmp(subjects_balloon{subject_balloon_id}, subjects_cardio));
   %-in 'combined_balloon' NaNs always appear together (column-wise), in 'combined_cardio' too
   if ~isempty(subject_cardio_id) && ~isnan(combined_balloon(2, 1, subject_balloon_id)) && ~isnan(combined_cardio(1, subject_cardio_id))
      %-left STG
      if ~isnan(combined_balloon(2, 1, subject_balloon_id)) && ~ismember(subject_balloon_id, outliers_balloon_aud)
         combined_balloon_cardio_aud(:, end+1) = [combined_balloon(2, 1:7, subject_balloon_id) combined_cardio(1:4, subject_cardio_id)'];
      end
   end
   if ~isempty(subject_cardio_id) && ~isnan(combined_balloon(4, 1, subject_balloon_id)) && ~isnan(combined_cardio(1, subject_cardio_id))
      %-left calcarine cortex
      if ~isnan(combined_balloon(4, 1, subject_balloon_id)) && ~ismember(subject_balloon_id, outliers_balloon_vis)
         combined_balloon_cardio_vis(:, end+1) = [combined_balloon(4, 1:7, subject_balloon_id) combined_cardio(1:4, subject_cardio_id)'];
      end
   end
end

%-auditory cortex: correlations between balloon and cardio
disp(['   '; 'aud'; '   ']);
for par_balloon_id     = 1:length(pars_balloon_t)
   for par_cardio_id   = 1:length(pars_cardio_t)
      [corr_R, corr_p] = corr(combined_balloon_cardio_aud(par_balloon_id, :)', combined_balloon_cardio_aud(7+par_cardio_id, :)');
      disp([pars_balloon_t{par_balloon_id} ' ' pars_cardio_t{par_cardio_id} ' ' num2str(corr_p)]);
   end
end

disp(['   '; 'vis'; '   ']);
%-visual cortex: correlations between balloon and cardio
for par_balloon_id     = 1:length(pars_balloon_t)
   for par_cardio_id   = 1:length(pars_cardio_t)
      [corr_R, corr_p] = corr(combined_balloon_cardio_vis(par_balloon_id, :)', combined_balloon_cardio_vis(7+par_cardio_id, :)');
      disp([pars_balloon_t{par_balloon_id} ' ' pars_cardio_t{par_cardio_id} ' ' num2str(corr_p)]);
   end
end

%-making figures
for cortex_id   =  1:2
   if cortex_id == 1
      combined_balloon_cardio = combined_balloon_cardio_aud;
      cortex = 'aud';
   else
      combined_balloon_cardio = combined_balloon_cardio_vis;
      cortex = 'vis';
   end
   figure('rend', 'painters', 'pos', [0 0 450 466], 'Visible', 'off');
   study_id_plot = 0;
   for par_balloon_id = 1:length(pars_balloon_t)
      for par_cardio_id  = 1:length(pars_cardio_t)
         study_id_plot = study_id_plot + 1;
         subplot(7, 4, study_id_plot);
         x_val = combined_balloon_cardio(7+par_cardio_id,   :)';
         y_val = combined_balloon_cardio(par_balloon_id, :)';
         axis off;
         %-[left bottom width height]
         gca_pos = get(gca, 'Position');
         %-lowering vertical spacing between subplots
         if par_balloon_id > 1
            gca_pos(2) = gca_pos(2) + (par_balloon_id-1)*0.004;
         end
         %-lowering horizontal spacing between subplots
         %-if sth changes in the figure, 0.043 might need to be adjusted...
         gca_pos(1)      = gca_pos(1) - (par_cardio_id-1)*0.043;
         axes('Position', gca_pos, 'Visible', 'off');
         par_for_subplot = combined_balloon_cardio(par_balloon_id, par_cardio_id, :);
         par_for_subplot = par_for_subplot(:);
         %-when using 'plot', circles are open...
         scatter(x_val, y_val, 2.5, 'filled', 'o'); hold on;
         %-without 'box on', only x and y axes are shown (no closed box)
         box on;
         axis([min(combined_balloon_cardio(7+par_cardio_id, :)) max(combined_balloon_cardio(7+par_cardio_id, :)) min(combined_balloon_cardio(par_balloon_id, :)) max(combined_balloon_cardio(par_balloon_id, :))]);
         lin_reg = fitlm(x_val, y_val);
         coefs   = lin_reg.Coefficients.Estimate;
         plot([min(x_val) max(x_val)], [coefs(1)+min(x_val)*coefs(2) coefs(1)+max(x_val)*coefs(2)], 'LineWidth', 1.25);
         [corr_R, corr_p] = corr(x_val, y_val);
         set(gca, 'FontSize', 4);
         if par_balloon_id < 7
            set(gca, 'XTickLabel', [], 'FontSize', 4);
         end
         if par_cardio_id > 1
            set(gca, 'YTickLabel', [], 'FontSize', 4);
         end
         title(['R=' sprintf('%.3f', corr_R) ', p=' sprintf('%.3f', corr_p)], 'FontSize', 5, 'FontWeight', 'normal');
         if par_cardio_id == 2
            text(1.024, 1.25, pars_balloon_t{par_balloon_id}, 'units', 'normalized', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
         end
         %-x-axis title: only for the last parameter -> last row of subplots
         if par_balloon_id == 7
            xlabel(pars_cardio_t{par_cardio_id}, 'FontSize', 5, 'FontWeight', 'bold');
         end
      end
   end
   figname = [paper '_balloon_cardio_' cortex];
   print_to_svg_to_pdf(figname, path_manage);
end
