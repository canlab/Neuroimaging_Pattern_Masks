clear all; close all;

addpath(genpath('~/software/spm12'));

addpath(genpath('~/software/canlab/CanlabCore'))
addpath(genpath('~/software/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('~/software/canlab/MasksPrivate'))

atlas_name = 'julich_fsl6';
space_description = 'MNI152NLin6Asym';
references = 'Amunts K, Mohlberg H, Blubau S, Zilles K. (2020) Julich_Brain: A 3D probablistic atlas of the human brain''s cytoarchitecture. Science 369(6506), 988-992.';

% imort atlas file in MNI152NLin2009cAsym space
% we just use this as a stand in template that we'll modify later, since
% this is the wrong space
MNI152NLin2009cAsym_bilat = which('JulichBrainAtlas_3.0_areas_MPM_b_N10_nlin2ICBM152asym2009c_public_11035603b4744231e17e87fd8ebcaf1a.nii.gz');

juStruct = parseXML(which('JulichBrainAtlas_3.0_areas_MPM_b_N10_nlin2ICBM152asym2009c_public_11035603b4744231e17e87fd8ebcaf1a.xml'));
labels = {};
for i = 2:2:length(juStruct.Children(2).Children)
    labels{end+1} = juStruct.Children(2).Children(i).Children.Data;
end

labels_R = cellfun(@(x1)(['R_',x1]),labels,'UniformOutput',false);
labels_L = cellfun(@(x1)(['L_',x1]),labels,'UniformOutput',false);

% combine data with labels
ref_file = which('fsl6_hcp_template.nii');
ref = fmri_data(ref_file);
juData = fmri_data(MNI152NLin2009cAsym_bilat).resample_space(ref,'nearest').remove_empty;
[~,~,juData.dat] = unique(juData.dat,'stable');
juAtlas = atlas(juData, ...
    'atlas_name', atlas_name,...
    'labels',[labels_L, labels_R], ...
    'labels_2',[labels_L, labels_R], ... % save the full labes here because we'll modify the labels field later
    'space_description', space_description, ...
    'references',references, 'noverbose');
%juAtlas = juAtlas.select_atlas_subset(find(~contains(juAtlas.labels,'GapMap')));
juAtlas = juAtlas.replace_empty();

% import probability maps
pmap = zeros(size(juAtlas.dat,1), length(juAtlas.labels));
parfor i = 1:length(juAtlas.labels)
    areaName = regexprep(strrep(regexprep(regexprep(juAtlas.labels{i},' \(.*',''),',.*',''),' ','-'),'[RL]_','');
    switch areaName
        case 'Medial-Accumbens'
            areaName = 'AcbM';
        case 'Lateral-Accumbens'
            areaName = 'AcbL';
        case 'Fundus-of-Caudate-Nucleus'
            areaName = 'FuCd';
        case 'Fundus-of-Putamen'
            areaName = 'FuP';
    end
    side = lower(regexprep(juAtlas.labels{i},'(^[RL]).*','$1'));
    file = dir(['probabilistic_maps_pmaps_157areas/', ...
        areaName, '/', areaName, '_pmap_', side, '_*nlin2ICBM152asym6*nii.gz']);
    try
        pdata = fmri_data([file.folder, '/', file.name]);    
        system(['rm -f ', file.folder, '/', strrep(file.name,'.gz','')]);
    catch
        error('failed to open %s', [file.folder, '/', file.name]);
    end

    % all this data is resampled, so let's get rid of small values which
    % are likely to extend outside the mas area
    pdata.dat(pdata.dat < 1e-4) = 0;

    % the gapmap probabilitys are bogus, so let's reset them. We pick 0.2
    % because anything below that gets thresholded later.
    if contains(juAtlas.labels{i},'GapMap')
        pdata.dat(pdata.dat > 0.2) = 0.2;
    end
    pmap(:,i) = pdata.dat;
end
delete(gcp('nocreate'))

juAtlas.probability_maps = sparse(pmap);

% now that we've extracted all the probability maps, let's rename the
% labels so that they can also be used as variables
juAtlas.labels = regexprep(regexprep(regexprep(juAtlas.labels,' \(.*',''),'[,.*]',''),'[- ]','_');

% Threshold at probability 0.2 or greater and k = 3 voxels or greater
juAtlas = threshold(juAtlas, .2, 'k', 3);

pureJuAtlas = juAtlas.select_atlas_subset(find(~contains(juAtlas.labels_2,'GapMap')));

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(pureJuAtlas, 'unique', 'overlay', which('fsl6_hcp_template.nii'));
figure;
% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(juAtlas);
pureR = atlas2region(pureJuAtlas)

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'full2', 'overlay', which('fsl6_hcp_template.nii'));
    brighten(.6)
    
    o2 = montage(pureR, o2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

    
end

%% save object

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save([pwd, '/' savename], 'juAtlas');
    
end

%% write - this writes only the label image

if dosave
    
    savename = sprintf('%s_atlas_regions.img', atlas_name);
    juAtlas.fullpath = fullfile(pwd, savename);
    write(juAtlas,'overwrite');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

for i = 1:length(r)
    
    eval([juAtlas.labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save([pwd, '/' savename], 'r', 'region_names', juAtlas.labels{:});

%%
if dosave
    
    figure; han = isosurface(pureJuAtlas);
    
    cellfun(@(x1)set(x1,'FaceAlpha', .5), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end
