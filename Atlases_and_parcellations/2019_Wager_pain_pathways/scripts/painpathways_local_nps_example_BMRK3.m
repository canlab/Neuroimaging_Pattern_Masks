% Load test dataset
test_dat = load_image_set('pain');  % bmrk3 data

sid = test_dat.additional_info.subject_id;
temp = test_dat.additional_info.temperatures;
pain = test_dat.Y;

% Load regions with local patterns stored in them
% Contains pain_regions_nps pain_regions_siips pain_regions_pdm1
load pain_pathways_region_obj_with_local_patterns % in Neuroimaging_Pattern_Masks

%% Local NPS patterns
% ---------------------------------------------------------------

% Dot product metric
[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_nps, test_dat);

% Cosine sim metric
% Note: For NPS, dot product metric is much higher than cosine sim metric.
% Normalizing data in local regions by overall intensity may help test how
% much the pattern (vs. overall average activity) contributes, but using
% cosine_similarity will not behave the same with local regions as it does
% with whole-brain patterns.
%[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_nps, test_dat, 'cosine_similarity');

% Correlation metric: Removes region averages
% Tests local pattern correlation, whereas dot product includes both
% pattern correlation and overall activation/scale of signal
[~, local_pattern_responses_corr] = extract_data(pain_regions_nps, test_dat, 'correlation');


% Plot correlation of pattern response with region averages
% Should be positive in regions with uniform positive weights, neg with
% negative weights. 
k = length(pain_regions_nps);
rr = [];
for i = 1:k
    rr(i) = corr(regions_with_testdata(i).dat, local_pattern_responses{i});
end

names = strrep({pain_regions_nps.shorttitle}, '_', ' ');

create_figure('Barplot_correlation_with_region_avgs: NPS');
bar(rr);
ylabel('Pearson''s r');
title('Correlation between local pattern response and region average');
set(gca, 'XTick', 1:length(rr), 'XLim', [0 length(rr)+1], 'XTickLabel', names, 'XTickLabelRotation', 45);
%saveas(gcf, fullfile('figures', 'r_local_pattern_region_average_nps.png'));

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
        
        % pattern correlation (demeaned)
        r_stim_corr(j, i) = corr(local_pattern_responses_corr{i}(wh), temp(wh));
        r_pain_corr(j, i) = corr(local_pattern_responses_corr{i}(wh), pain(wh));
        
        % average activity
        r_stim_avg(j, i) = corr(regions_with_testdata(i).dat(wh), temp(wh));
        r_pain_avg(j, i) = corr(regions_with_testdata(i).dat(wh), pain(wh));
        
        
        
    end
end

% Line plots for correlations with temperature and pain
% -----------------------------------------------------------------------

color1 = [.7 .4 0];
color2 = [.2 .4 .7];
color3 = [.2 .7 .2];

create_figure('Correlations', 2, 1);

han = barplot_columns(r_stim_corr, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color1, 'MarkerFaceColor',color1, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color1, 'LineStyle', 'none');
set(han.star_handles, 'Color', color1);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .92; set(han.star_handles(i), 'Position', mypos); end

xlabel(' '); title('Avg within-person correlation with temperature'); 
ylabel('Within-person r');
hh = plot_horizontal_line(0); set(hh, 'LineStyle', '--');
set(gca, 'YLim', [-.5 1]);

lh1 = han.line_han(1);      % for legend

han = barplot_columns(r_stim_dot, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color2, 'MarkerFaceColor',color2, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color2, 'LineStyle', 'none');
set(han.star_handles, 'Color', color2);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .98; set(han.star_handles(i), 'Position', mypos); end

lh2 = han.line_han(1);

han = barplot_columns(r_stim_avg, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color3, 'MarkerFaceColor',color2, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color3, 'LineStyle', 'none');
set(han.star_handles, 'Color', color3);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .86; set(han.star_handles(i), 'Position', mypos); end

lh3 = han.line_han(1);

legend([lh2, lh1, lh3], {'Dot product' 'Pattern Correlation' 'Region average'});

subplot(2, 1, 2);
han = barplot_columns(r_pain_corr, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color1, 'MarkerFaceColor',color1, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color1, 'LineStyle', 'none');
set(han.star_handles, 'Color', color1);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .92; set(han.star_handles(i), 'Position', mypos); end

xlabel(' '); title('Avg within-person correlation with pain');
ylabel('Within-person r');
hh = plot_horizontal_line(0); set(hh, 'LineStyle', '--');
set(gca, 'YLim', [-.5 1]);

han = barplot_columns(r_pain_dot, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color2, 'MarkerFaceColor',color2, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color2, 'LineStyle', 'none');
set(han.star_handles, 'Color', color2);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .98; set(han.star_handles(i), 'Position', mypos); end

han = barplot_columns(r_stim_avg, 'names', names, 'nofigure', 'noviolin', 'noind', 'line');
set(han.line_han, 'Color', color3, 'MarkerFaceColor',color2, 'LineStyle', 'none');
set(han.errorbar_han, 'Color', color3, 'LineStyle', 'none');
set(han.star_handles, 'Color', color3);
for i = 1:length(han.star_handles), mypos = get(han.star_handles(i), 'Position'); mypos(2) = .86; set(han.star_handles(i), 'Position', mypos); end

set(gca, 'YLim', [-.5 1]);
