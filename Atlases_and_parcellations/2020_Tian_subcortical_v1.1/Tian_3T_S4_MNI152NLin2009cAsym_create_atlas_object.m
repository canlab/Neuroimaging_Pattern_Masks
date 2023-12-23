% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI
addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))

atlas_name = 'tian_3t_fsl6';
space_description = 'MNI152NLin6Asym';
references = 'Tian Y, Margulies D, Breakspear M, Zalesky A (2020). Nature Neuroscience. 23(11) 1421-1432.';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% Authors only provide 1mm parcellations in this space together with
% cortical parcellations, so we'll need to take this and truncate it.
parcellation_file = which('Schaefer2018_100Parcels_7Networks_order_Tian_Subcortex_S4_3T_MNI152NLin2009cAsym_1mm.nii.gz');


% Get labels
% -----------------------------------------------------------------------
fid = fopen(which('Tian_Subcortex_S4_3T_label.txt'));
labels = textscan(fid,'%s');
labels = cellfun(@(x1)strrep(x1,'-','_'),labels{1},'UniformOutput',false);
labels = labels(:)';
fclose(fid);

%% identify iid individuals
% this section doesn't do anything. It was just used to generate a list
% which was then saved and used on an HPC CLI to extract the necessary
% files from Tian2020MSAProbmap's.

% this list gives us the order of participants in the segmentation files
fid = fopen(which('Tian2020MSAProbmap/subject_test_list.txt'));
Tian_sid = textscan(fid,'%s');
fclose(fid);
Tian_sid  = Tian_sid{1};

% HCP data has a bunch of individuals from the same families. For unbiased
% probability maps we need the to sample independent individuals who are
% unrelated. This list contains those. Unfortuntaely you need permission
% from HCP to access these lists, and while I have it, I don't have
% permission to share these files. You'll have to obtain access to the
% restrictied subject demographic data from HCP for yourself to get these
% files.

tbl = readtable(['/home/bogdan/MyDocuments/canlab/hcp_hyperalignment/', '/docs/RESTRICTED_bogpetre_9_12_2023_9_35_28.csv']);
utbl = readtable(['/home/bogdan/MyDocuments/canlab/hcp_hyperalignment/', '/docs/unrestricted_bogpetre_9_14_2023_8_25_4.csv']);

% let's filter tables for the subject Tian used
tbl = tbl(ismember(tbl.Subject,cellfun(@str2num,Tian_sid)),:);
utbl = utbl(ismember(tbl.Subject,cellfun(@str2num,Tian_sid)),:);

% goal is to maximize number of independent families, then maximize number
% of twins over unrelated siblings
blacklist = ismember(tbl.Subject,[143527, ... % This subject shares a mother with one individual and a father with another individual who are otherwise independent of one another
    211316, ... % is half sibling with a DZ pair on mother side and half sibling with an independent individual on fathers side who is otherwise independent
    713239]); % is half sibling with a DZ pair on mother side and half sibling with an independent individual on fathers side who is otherwise independent


uniq_fid = unique(tbl.Family_ID(~blacklist));

% pick an individual from each family. Prefer genotypically verified twins
iid_individ = zeros(height(tbl),1);
for i = 1:length(uniq_fid)
    this_tbl = tbl(ismember(tbl.Family_ID, uniq_fid(i)),:);
    if any(contains(this_tbl.ZygosityGT,{'MZ','DZ'}))
        this_tbl = this_tbl(contains(this_tbl.ZygosityGT,{'MZ','DZ'}),:);
    end
    this_sid = this_tbl.Subject(randperm(height(this_tbl),1));
    iid_individ = iid_individ | tbl.Subject == this_sid;
end

iid_individ_ind = find(iid_individ);

keep_sid = tbl.Subject(iid_individ_ind);
keep_ind = find(ismember(cellfun(@str2num,Tian_sid),keep_sid));

fid=fopen('/tmp/iid_sids.csv','w+'); 
fprintf(fid, '%d\n', keep_ind); 
fclose(fid)

%% import individual subject segmentations
% truncate cortical areas
parcellation = fmri_data(parcellation_file);
parcellation.dat(parcellation.dat > length(labels)) = 0;

fnames = dir('iid_parcellations/*nii.gz');
pmaps = cell(length(fnames),1);
parfor (f = 1:length(fnames),5)
    this_pmaps = fmri_data(fullfile(fnames(f).folder, fnames(f).name));
    pmaps{f} = atlas(this_pmaps);
end

fmaps = cellfun(@fmri_data, pmaps, 'UniformOutput', false);
all_maps = cat(fmaps{:});

parcels = zeros(size(all_maps.dat,1),num_regions(pmaps{1}));
for i = 1:num_regions(pmaps{1})
    this_region = (all_maps.dat == i);
    parcels(:,i) = mean(this_region,2);
end
new_map = all_maps.get_wh_image(1); 
new_map.dat = parcels;

new_atlas = atlas(new_map, 'labels', labels).threshold(0,'k',5);


new_atlas = atlas(new_map, 'labels', labels);

% the above values aren't probabilities. They don't sum to 1.
% we can sort of assume they sum to 1 within the relevant CIFTI mask, 
% since they were generated from HCP CIFTI files these are from the 
% 2023_CANLab_atlas/src folder
cifti_subctx = fmri_data(which('hcp_cifti_subctx_labels.nii.gz'));
fid = fopen(which('hcp_cifti_subctx_labels.txt'),'r');
cifti_lbls = textscan(fid,'%s\n');
fclose(fid);
cifti_subctx_atlas = atlas(cifti_subctx, 'labels', cifti_lbls{1}');
tian_regions = {'accumbens', 'amygdala', 'caudate', 'hippocampus', 'pallidum', 'putamen', 'thalamus'};
cifti_mask = fmri_mask_image(cifti_subctx_atlas.select_atlas_subset(tian_regions));
cifti_mask = cifti_mask.resample_space(pmaps,'nearest');
pmaps = pmaps.apply_mask(cifti_mask).remove_empty();




pmaps.fullpath = 'Tian_3T_S4_MNI152NLin6Asym_probability_maps.nii';
pmaps.write();
gzip(pmaps.fullpath)
delete(pmaps.fullpath)

% Create object
% -----------------------------------------------------------------------
atlas_obj = atlas(pmaps, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'space_description', space_description, ...
    'references', references, 'noverbose');

atlas_obj.threshold(0,'k',5,'remove_parcel_fragments').orthviews

% I've emailed Tian and Zalesky asking about probability maps, they have
% them in figure 7c of their paper, but for now we default to a binary map.
uniq_id = unique(atlas_obj.remove_empty.dat);
n_roi = length(uniq_id);
pmap = zeros(atlas_obj.volInfo.n_inmask, n_roi);
atlas_obj = atlas_obj.replace_empty();
for i = 1:length(uniq_id)
    pmap(:,i) = atlas_obj.dat == uniq_id(i);
end
atlas_obj.probability_maps = pmap;


% Process object
% -----------------------------------------------------------------------

% Threshold at probability 0.2 or greater and k = 3 voxels or greateratlas_obj = threshold(atlas_obj, 0.2, 'k', 3);
atlas_obj = threshold(atlas_obj, .2, 'k', 3);


% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
atlas_obj.orthviews('unique',which('MNI152NLin2009cAsym_T1_1mm.nii.gz'));


% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'multirow', 1, 'overlay', which('MNI152NLin2009cAsym_T1_1mm.nii.gz'));
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

    
end
 
%% save object

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save(savename, 'atlas_obj');
    
end

%% write - this writes only the label image

if dosave
    
    savename = sprintf('%s_atlas_regions.img', atlas_name);
    atlas_obj.fullpath = fullfile(pwd, savename);
    write(atlas_obj,'overwrite');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    cellfun(@(x1)(set(x1,'FaceAlpha', .5)), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end
