addpath(genpath('~/.matlab/canlab/Neuroimaging_Pattern_Masks'));
addpath(genpath('~/.matlab/canlab/CanlabCore'));
addpath('~/.matlab/spm/spm12');

% load data
indAtlPath = '/home/bogdan/MyDocuments/canlab/atlases/results/thalamus/spacetop';

% the following are "simple" contrasts from heejungs first level designs
% for spacetop. The files were obtained from discovery by concatenating
% contrasts from 
% /dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_atlas/stats/spacetop
% for subjects listed in
% Atlases_and_parcellations/2018_Iglesias_thalamic_HCP278_reconstruction/src/spacetop_participants.csv
pain = fmri_data('/home/bogdan/MyDocuments/canlab/atlases/sandbox/thalamic_seg_demo/con_0013.nii.gz');
vic = fmri_data('/home/bogdan/MyDocuments/canlab/atlases/sandbox/thalamic_seg_demo/con_0014.nii.gz');
cog = fmri_data('/home/bogdan/MyDocuments/canlab/atlases/sandbox/thalamic_seg_demo/con_0015.nii.gz');
sid = readtable(which('spacetop_participants.csv'),'ReadVariableNames',false);

thal_atlas = load_atlas('iglesias_thal_fmriprep20').threshold(0.2);

% replace thalamic atlas probabilities with spacetop specific ones
spacetop76 = fmri_data(which('spacetop76b.ThalamicNuclei.v13.T1DWI.MNI152NLin2009cAsym.nii.gz'));
has_contrast = ~ismember(sid.Var1,{'sub-0017','sub-0063','sub-0075'});
spacetop73 = spacetop76.get_wh_image(find(has_contrast));
sid = sid(has_contrast,:);
    
n_regions = length(unique(spacetop73.dat)) - 1;
uniq_rois = unique(spacetop73.dat);
uniq_rois(uniq_rois == 0) = [];

pmap = zeros(size(spacetop73.dat,1), n_regions);
for i = 1:n_regions
    this_roi = uniq_rois(i);
    this_mask = (spacetop73.dat == this_roi);
    pmap(:,i) = mean(this_mask,2);
end

thal_atlas.probability_maps = pmap;
thal_atlas = thal_atlas.probability_maps_to_region_index;

% clean up the new atlas to get rid of some stray labels in brainstem
mask = fmri_mask_image(thal_atlas.threshold(3/76)).replace_empty();
mask.dat = iimg_smooth_3d(mask.dat,mask.volInfo,3);
thal_atlas = thal_atlas.apply_mask(mask);

% threshold for a more representative atlas
thal_atlas = thal_atlas.threshold(0.2);

% extract mean ROI contrast based on group atlas
thal_atlas = thal_atlas.resample_space(pain);
pain_mean_grpAtl = thal_atlas.extract_data(pain);
vic_mean_grpAtl = thal_atlas.extract_data(vic);
cog_mean_grpAtl = thal_atlas.extract_data(cog);

% extract mean ROI contrast based on individual atlases
indAtl = spacetop73;

fsind = [8103:8106,8108:8113,8115,8116,8118,8120:8122,8126:8129,8133,8135,8136,8203:8206,8208:8213,8215,8216,8218,8220:8222,8226:8229,8233,8235,8236];
[pain_mean_indAtl, vic_mean_indAtl, cog_mean_indAtl] = deal(nan(height(sid),length(thal_atlas.labels)));

parfor i = 1:size(spacetop73.dat,2)
    this_atlas = indAtl.get_wh_image(i);
    this_atlas = this_atlas.resample_space(pain.get_wh_image(i),'nearest');
    
    % identify which regions are present (some may be missing)
    this_atlas_regions = 1:46;
    for ind = 1:length(fsind)
        if length(this_atlas.dat(this_atlas.dat == fsind(ind))) == 0
            this_atlas_regions(ind) = [];
        end
    end
    this_atlas = atlas(this_atlas,...
        'labels', thal_atlas.labels(this_atlas_regions));
    this_atlas = this_atlas.apply_mask(mask);
    this_atlas = this_atlas.check_properties('compress_index');

    % extract data and assign to indices corresponding to regions we have
    [this_pain_mean_indAtl, this_vic_mean_indAtl, this_cog_mean_indAtl] = deal(nan(1,length(thal_atlas.labels)));
    this_pain_mean_indAtl(this_atlas_regions) = this_atlas.extract_data(pain.get_wh_image(i));
    this_vic_mean_indAtl(this_atlas_regions) = this_atlas.extract_data(vic.get_wh_image(i));
    this_cog_mean_indAtl(this_atlas_regions) = this_atlas.extract_data(cog.get_wh_image(i));

    % the above operation can't be performed directly on output vars in a
    % parfor loop, so now we assign results to our put vars.
    pain_mean_indAtl(i,:) = this_pain_mean_indAtl;
    vic_mean_indAtl(i,:) = this_vic_mean_indAtl;
    cog_mean_indAtl(i,:) = this_cog_mean_indAtl;
end

%% stats

% make long form table
% stack regions
% subjects nested within regions
subjects = repmat((1:size(pain_mean_grpAtl,1))',1,size(pain_mean_grpAtl,2));
% dummy code regions
regions = {};
for i = 1:size(pain_mean_grpAtl,2)
    regions{i} = zeros(size(pain_mean_grpAtl,1), size(pain_mean_grpAtl,2));
    regions{i}(:,i) = 1;
end

pain_grpAtl_vert = reshape(pain_mean_grpAtl,prod(size(pain_mean_grpAtl)),1);
vic_grpAtl_vert = reshape(vic_mean_grpAtl,prod(size(vic_mean_grpAtl)),1);
cog_grpAtl_vert = reshape(cog_mean_grpAtl,prod(size(cog_mean_grpAtl)),1);

pain_indAtl_vert = reshape(pain_mean_indAtl,prod(size(pain_mean_indAtl)),1);
vic_indAtl_vert = reshape(vic_mean_indAtl,prod(size(vic_mean_indAtl)),1);
cog_indAtl_vert = reshape(cog_mean_indAtl,prod(size(cog_mean_indAtl)),1);

regions_vert = cellfun(@(r)reshape(r,prod(size(r)),1), regions, 'UniformOutput',false);
subjects_vert = reshape(subjects,prod(size(subjects)),1);

% stack tasks
% subjects nested within regions nested within tasks
indAtl_vert = [pain_indAtl_vert; vic_indAtl_vert; cog_indAtl_vert];
grpAtl_vert = [pain_grpAtl_vert; vic_grpAtl_vert; cog_grpAtl_vert];
pain_v_vic = reshape(repmat([0.5,-0.5,0],size(pain_indAtl_vert,1),1),3*size(pain_indAtl_vert,1),1);
pain_vic_v_cog = reshape(repmat([-1/3,-1/3,2/3],size(pain_indAtl_vert,1),1),3*size(pain_indAtl_vert,1),1);
regions_vert3 = cellfun(@(r)repmat(r,3,1), regions_vert, 'UniformOutput',false);
subjects_vert3 = repmat(subjects_vert,3,1);

% stack individ and group parcellations
% subjects nested within regions nested within tasks
individ_parcellation = 0.5*ones(size(indAtl_vert));
group_parcellation = -0.5*ones(size(grpAtl_vert));

dataArr = [double([indAtl_vert; grpAtl_vert]),...
    double([pain_v_vic; pain_v_vic]),...
    double([pain_vic_v_cog; pain_vic_v_cog]),...
    [cat(2,regions_vert3{:}); cat(2,regions_vert3{:})],...
    double([subjects_vert3; subjects_vert3]),...
    double([individ_parcellation; group_parcellation])];
tbl = array2table(dataArr, ...
    'VariableNames', ...
    ['BOLD','pain_v_vic','pain_vic_v_cog',thal_atlas.labels,'subject','parcellation']);

[stats, mdl] = deal(cell(num_regions(thal_atlas),1));
parfor i = 1:num_regions(thal_atlas)
    [stats{i}, mdl{i}] = fit_single_region_model(tbl,thal_atlas.labels(i));
end

for i = 1:length(stats)
    disp(thal_atlas.label_descriptions{i});
    disp(stats{i});
    disp('-----------------------------');
end

%% post-hoc tests of main effects
% The following plots show regions where the main effect of region was
% significantly different between individual and group parcellation. The
% three regions where contrast was most enhanced and most decreased are
% shown.

[~,I] = sort(cellfun(@(x1)(x1.tStat(4)), stats));
I = [I(1:3); I(end-2:end)];
for i = 1:length(I)
    this_ind = I(i);
    plot_post_hoc_contrast(tbl,thal_atlas,thal_atlas.labels{this_ind}, 0.9);
    sgtitle(sprintf('%s (t=%0.2f, p = %0.3e)', strrep(thal_atlas.label_descriptions{this_ind},'_',' '),...
        stats{this_ind}.tStat(4), stats{this_ind}.pValue(4)));
end

%% post-hoc tests of task contrasts
% the following plots show regions where between task contrasts
% significantly differ between group and individual parcellations
% the two regions regions where between task contrasts are most improved
% and most degraded are plotted.

[~,I] = sort([cellfun(@(x1)(x1.tStat(5)), stats); cellfun(@(x1)(x1.tStat(6)), stats)]);
I = [I(1:2); I(end-1:end)];
for i = 1:length(I)
    this_roi_ind = I(i);
    if this_roi_ind <= length(stats)
        t = stats{this_roi_ind}.tStat(5);
        p = stats{this_roi_ind}.pValue(5);
        con = 'pain vs. vic';
    else
        this_roi_ind = this_roi_ind - length(stats);
        t = stats{this_roi_ind}.tStat(6);
        p = stats{this_roi_ind}.pValue(6);
        con = 'pain & vic vs cog';
    end
    plot_post_hoc_contrast(tbl,thal_atlas,thal_atlas.labels{this_ind}, 0.9);
    sgtitle(sprintf('%s (%s: t=%0.2f, p = %0.3e)', strrep(thal_atlas.label_descriptions{this_roi_ind},'_',' '),con,t,p));
end

function fig = plot_post_hoc_contrast(tbl, atlas_obj, roi_label, border_prop)

tbl_roi = tbl(tbl.(roi_label) == 1,:);
ind_roi_pain = tbl_roi.BOLD(tbl_roi.pain_v_vic == 0.5 & tbl_roi.parcellation > 0);
ind_roi_vic = tbl_roi.BOLD(tbl_roi.pain_v_vic == -0.5 & tbl_roi.parcellation > 0);
ind_roi_cog = tbl_roi.BOLD(tbl_roi.pain_vic_v_cog > 0 & tbl_roi.parcellation > 0);

grp_roi_pain = tbl_roi.BOLD(tbl_roi.pain_v_vic == 0.5 & tbl_roi.parcellation < 0);
grp_roi_vic = tbl_roi.BOLD(tbl_roi.pain_v_vic == -0.5 & tbl_roi.parcellation < 0);
grp_roi_cog = tbl_roi.BOLD(tbl_roi.pain_vic_v_cog > 0 & tbl_roi.parcellation < 0);

roi = atlas2region(atlas_obj.threshold(0.2).select_atlas_subset({roi_label},'exact'));
[f,ax] = deal(cell(3,1));
orientations = {'saggital','axial','coronal'};
for i = 1:length(orientations)
    ax{i} = roi.montage('regioncenters','nofigure',orientations{i});
    f{i} = gcf;
    set(f{i},'Tag',orientations{i})
end

fig = figure;
t0 = tiledlayout(1,3);
t1 = tiledlayout(t0,2,2,'Padding','compact','TileSpacing','compact');
ax1 = nexttile(t1);
copyobj(ax{1}.montage{1}.axis_handles.Children, ax1);
colormap('gray');
axis square
axis off;

xyzminmax = [min(roi.XYZmm, [], 2) max(roi.XYZmm, [], 2)];  % 3 x 2, 3 min  3 max mm coords
xl = xyzminmax(2, :);
yl = xyzminmax(3, :);
myborder = border_prop * range(xl);     % multiple of range from min to max coordinate
myborder = max([myborder 20]);          % border at least 10 mm
xl(1) = xl(1) - myborder;
xl(2) = xl(2) + myborder;
myborder = border_prop * range(yl);
myborder = max([myborder 20]);          % border at least 10 mm
yl(1) = yl(1) - myborder;
yl(2) = yl(2) + myborder;
xlim(xl);
ylim(yl);
close(f{1})

ax2 = nexttile(t1);
copyobj(ax{2}.montage{1}.axis_handles.Children, ax2);
colormap('gray');
axis square
axis off;

xyzminmax = [min(roi.XYZmm, [], 2) max(roi.XYZmm, [], 2)];  % 3 x 2, 3 min  3 max mm coords
xl = xyzminmax(1, :);
yl = xyzminmax(2, :);
myborder = border_prop * range(xl);     % multiple of range from min to max coordinate
myborder = max([myborder 20]);          % border at least 10 mm
xl(1) = xl(1) - myborder;
xl(2) = xl(2) + myborder;
myborder = border_prop * range(yl);
myborder = max([myborder 20]);          % border at least 10 mm
yl(1) = yl(1) - myborder;
yl(2) = yl(2) + myborder;
xlim(xl);
ylim(yl);

close(f{2});

ax3 = nexttile(t1); axis off
ax4 = nexttile(t1)
copyobj(ax{3}.montage{1}.axis_handles.Children, ax4);
colormap('gray');
axis square
axis off;

xyzminmax = [min(roi.XYZmm, [], 2) max(roi.XYZmm, [], 2)];  % 3 x 2, 3 min  3 max mm coords
xl = xyzminmax(1, :);
yl = xyzminmax(3, :);
myborder = border_prop * range(xl);     % multiple of range from min to max coordinate
myborder = max([myborder 20]);          % border at least 10 mm
xl(1) = xl(1) - myborder;
xl(2) = xl(2) + myborder;
myborder = border_prop * range(yl);
myborder = max([myborder 20]);          % border at least 10 mm
yl(1) = yl(1) - myborder;
yl(2) = yl(2) + myborder;
xlim(xl);
ylim(yl);
close(f{3})

axM = nexttile(t0)
barplot_columns([ind_roi_pain, ind_roi_vic, ind_roi_cog],'nofig');
title('Individualized Atlas','FontSize',14);
ylabel('ROI mean \beta');
xlabel('task');
set(gca,'XTickLabels',{'Pain','Vic','Cog'})
grid on;

axR = nexttile(t0)
barplot_columns([grp_roi_pain, grp_roi_vic, grp_roi_cog],'nofig');
title('Group Mean (N=73) Atlas','FontSize',14)
ylabel('ROI mean \beta');
xlabel('task');
set(gca,'XTickLabels',{'Pain','Vic','Cog'})
grid on;

yl1 = get(axM,'YLim');
yl2 = get(axR,'YLim');
yl = [min(yl1(1),yl2(1)), max(yl1(2),yl2(2))];
set(axM,'YLim',yl);
set(axR,'YLim',yl);

pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),1000,480])
sgtitle(strrep(roi_label,'_',' '));

end

function [stats, mdl] = fit_single_region_model(tbl,region_labels)
    [mdl_str1, mdl_str2, mdl_str3] = deal([]);
    for i = 1:length(region_labels)
        mdl_str1 = [mdl_str1, sprintf(' + %s',region_labels{i})];
        mdl_str2 = [mdl_str2, sprintf(' + pain_v_vic:%s',region_labels{i})];
        mdl_str3 = [mdl_str3, sprintf(' + pain_vic_v_cog:%s',region_labels{i})];
    end

    [mdl_str4, mdl_str5, mdl_str6] = deal([]);
    for i = 1:length(region_labels)
        mdl_str4 = [mdl_str4, sprintf(' + parcellation:%s',region_labels{i})];
        mdl_str5 = [mdl_str5, sprintf(' + parcellation:pain_v_vic:%s',region_labels{i})];
        mdl_str6 = [mdl_str6, sprintf(' + parcellation:pain_vic_v_cog:%s',region_labels{i})];
    end
    mdl_str = [mdl_str1, mdl_str2, mdl_str3, mdl_str4, mdl_str5, mdl_str6];
    mdl_str = mdl_str(3:end);

    new_tbl = cell(1,length(region_labels));
    for i = 1:length(region_labels)
        new_tbl{i} = tbl(tbl.(region_labels{i}) == 1,:);
    end
    tbl = cat(1,new_tbl{:});

    mdl = fitlme(tbl,sprintf('BOLD ~ %s -1 + (%s -1 | subject)',mdl_str,mdl_str),'CovariancePattern','full','FitMethod','REML');
    [~,~,stats] = fixedEffects(mdl,'dfmethod','satterthwaite');
end