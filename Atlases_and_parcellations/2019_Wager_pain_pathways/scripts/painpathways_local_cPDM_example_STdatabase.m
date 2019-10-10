
cd('/Users/torwager/Documents/GitHub/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2019_Wager_pain_pathways')

load data/pain_pathways_single_trial_data_2019 ST_cleaned
load pain_pathways_atlas_obj.mat

% Load regions with local patterns stored in them
% Contains pain_regions_nps pain_regions_siips pain_regions_pdm1
load pain_pathways_region_obj_with_local_patterns % in Neuroimaging_Pattern_Masks

sid = ST_cleaned.subj_idx;
temp = ST_cleaned.abs_temp;
pain = ST_cleaned.pain_rating;


%% Local cPDM patterns
% ---------------------------------------------------------------

% Dot product metric - already loaded
%[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_nps, test_dat);
local_pattern_responses = {};

k = size(ST_cleaned.cPDM, 2);

for i = 1:k
    
    local_pattern_responses{i} = ST_cleaned.cPDM(:, i);
    
end

names = strrep(ST_cleaned.cpdm_labels, '_', ' ');

% % regions_with_testdata = ****

% % Plot correlation of pattern response with region averages
% % Should be positive in regions with uniform positive weights, neg with
% % negative weights. 
% k = length(pain_regions_nps);
% rr = [];
% for i = 1:k
%     rr(i) = corr(regions_with_testdata(i).dat, local_pattern_responses{i});
% end


% % create_figure('Barplot_correlation_with_region_avgs: NPS');
% % bar(rr);
% % ylabel('Pearson''s r');
% % title('Correlation between local pattern response and region average');
% % set(gca, 'XTick', 1:length(rr), 'XLim', [0 length(rr)+1], 'XTickLabel', names, 'XTickLabelRotation', 45);
% % %saveas(gcf, fullfile('figures', 'r_local_pattern_region_average_nps.png'));

% Average (across subjects) within-person correlation with stim intensity and pain
n = max(sid);

% Calculate correlations for both the dotproduct and local Pearson's r metrics
% -----------------------------------------------------------------------

[r_stim_dot, r_pain_dot, r_stim_corr, r_pain_corr] = deal(zeros(n, k));

for i = 1:k
    for j = 1:n
        
        wh = sid == j;
        
        % dot product
        r_stim_dot(j, i) = corr(local_pattern_responses{i}(wh), temp(wh));
        r_pain_dot(j, i) = corr(local_pattern_responses{i}(wh), pain(wh));
        
%         % pattern correlation (demeaned)
%         r_stim_corr(j, i) = corr(local_pattern_responses_corr{i}(wh), temp(wh));
%         r_pain_corr(j, i) = corr(local_pattern_responses_corr{i}(wh), pain(wh));
%         
%         % average activity
%         r_stim_avg(j, i) = corr(regions_with_testdata(i).dat(wh), temp(wh));
%         r_pain_avg(j, i) = corr(regions_with_testdata(i).dat(wh), pain(wh));
%         
        
        
    end
end

% Line plots for correlations with temperature and pain
% -----------------------------------------------------------------------

color1 = [.7 .4 0];
color2 = [.2 .4 .7];
color3 = [.2 .7 .2];

create_figure('Correlations', 2, 1);

% han = barplot_columns(r_stim_corr, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
% set(han.line_han, 'Color', color1, 'MarkerFaceColor',color1, 'LineStyle', 'none');
% set(han.errorbar_han, 'Color', color1, 'LineStyle', 'none');
% set(han.star_handles, 'Color', color1);
% for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .92; set(han.star_handles(i), 'Position', mypos); end

% lh1 = han.line_han(1);      % for legend

han = barplot_columns(r_stim_dot, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color2, 'MarkerFaceColor',color2, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color2, 'LineStyle', 'none');
set(han.star_handles, 'Color', color2);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .98; set(han.star_handles(i), 'Position', mypos); end

hh = plot_horizontal_line(0); set(hh, 'LineStyle', '--'); 
xlabel(' '); title('Avg within-person correlation (single trials)'); 
ylabel('Single-trial r (temperature)');
set(gca, 'YLim', [-0.05 .2]);

lh2 = han.line_han(1);

% han = barplot_columns(r_stim_avg, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
% set(han.line_han, 'Color', color3, 'MarkerFaceColor',color2, 'LineStyle', 'none');
% set(han.errorbar_han, 'Color', color3, 'LineStyle', 'none');
% set(han.star_handles, 'Color', color3);
% for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .86; set(han.star_handles(i), 'Position', mypos); end

% lh3 = han.line_han(1);

% legend([lh2, lh1, lh3], {'Dot product' 'Pattern Correlation' 'Region average'});

subplot(2, 1, 2);
% han = barplot_columns(r_pain_corr, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
% set(han.line_han, 'Color', color1, 'MarkerFaceColor',color1, 'LineStyle', 'none');
% set(han.errorbar_han, 'Color', color1, 'LineStyle', 'none');
% set(han.star_handles, 'Color', color1);
% for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .92; set(han.star_handles(i), 'Position', mypos); end


han = barplot_columns(r_pain_dot, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color2, 'MarkerFaceColor',color2, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color2, 'LineStyle', 'none');
set(han.star_handles, 'Color', color2);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .98; set(han.star_handles(i), 'Position', mypos); end

% han = barplot_columns(r_pain_avg, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
% set(han.line_han, 'Color', color3, 'MarkerFaceColor',color2, 'LineStyle', 'none');
% set(han.errorbar_han, 'Color', color3, 'LineStyle', 'none');
% set(han.star_handles, 'Color', color3);
% for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .86; set(han.star_handles(i), 'Position', mypos); end

xlabel(' '); %title('Avg within-person correlation with pain');
ylabel('Single-trial r (pain)');
set(gca, 'YLim', [-0.05 .2]);
hh = plot_horizontal_line(0); set(hh, 'LineStyle', '--'); 

saveas(gcf, 'painpathways_cPDM_temp_pain_relationships.png');

