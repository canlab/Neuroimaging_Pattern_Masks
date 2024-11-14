clear all;
close all;

addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('/home/bogdan/.matlab/canlab/MasksPrivate'))

leaddbsRoot='/home/bogdan/MyDocuments/canlab/atlases/_resources/LeadsDBS/';

myAtlas = load_atlas('iglesias_thal_fmriprep20');

% note, this requires MasksPrivate
morel = load_atlas('morel_fmriprep20');

%% plot overlaps

% equivalence mapping {{iglesias}, {morel}}
emap = {{{'AV'},{'AD','AM','AV'}},...
    {{'CeM'},{'CeM'}},...
    {{'CL'},{'CL'}},...
    {{'CM'},{'CM'}},...
    {{'LD'},{'LD'}},...
    {{'LGN'},{'LGNmc','LGNpc'}},...
    {{'LP'},{'LP'}},...
    {{'L_SG'},{'SG','Li'}},...
    {{'MDl'},{'MDpc'}},...
    {{'MDm'},{'MDmc'}},...
    {{'MGN'},{'MGN'}},...
    {{'MV_Re'},{'MV'}},...
    {{'Pf'},{'Pf'}},...
    {{'PuA'},{'PuA'}},...
    {{'PuI'},{'PuI','Po'}},...
    {{'PuL'},{'PuL'}},...
    {{'VA'},{'VApc'}},...
    {{'VAmc'},{'VAmc'}},...
    {{'VLa'},{'VLa'}},...
    {{'VLp'},{'VLpd','VLpv','VM'}},...
    {{'VPL_VPM'},{'VPLa','VPLp','VPI','VPM'}},...
    {{'PuMm','PuMl'},{'PuM'}}};

thisAtlas = myAtlas.select_atlas_subset({'L_'}).threshold(0.2);
thisAtlas = thisAtlas.select_atlas_subset(find(~contains(thisAtlas.labels,{'R_L'})));
thisAtlas.probability_maps = [];
thisMorelAtlas = morel.select_atlas_subset({'L_'});

[newAtlas1, newAtlas2] = deal({});
for i = 1:length(emap)
    newAtlas1{end+1} = thisAtlas.select_atlas_subset(emap{i}{1},'flatten');
    newAtlas2{end+1} = thisMorelAtlas.select_atlas_subset(emap{i}{2},'flatten');
end
thisAtlas = [newAtlas1{:}];
thisMorelAtlas = [newAtlas2{:}];

for orientation = {'saggital','coronal','axial'}
    %%
    o2 = thisAtlas.montage('transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(thisAtlas)
        try
            morel_roi = thisMorelAtlas.select_atlas_subset(i);
        
            if num_regions(morel_roi) == 1
                o3 = o2;
                o3.activation_maps = o2.activation_maps(i);
                o3.montage = o2.montage(i);
                morel_roi.montage(o3,'existing_figure','existing_axes', o2.montage{i}.axis_handles,'outline','color',[0,0,0]);
            end
        end
    end
    
    set(gcf,'Tag',orientation{1});
end


thisAtlas = myAtlas.select_atlas_subset({'R_'}).threshold(0.2);
thisAtlas.probability_maps = [];
thisMorelAtlas = morel.select_atlas_subset({'R_'});

[newAtlas1, newAtlas2] = deal({});
for i = 1:length(emap)
    newAtlas1{end+1} = thisAtlas.select_atlas_subset(emap{i}{1},'flatten');
    newAtlas2{end+1} = thisMorelAtlas.select_atlas_subset(emap{i}{2},'flatten');
end
thisAtlas = [newAtlas1{:}];
thisMorelAtlas = [newAtlas2{:}];

for orientation = {'saggital','coronal','axial'}
    %%
    o2 = thisAtlas.montage('transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(thisAtlas)
        try
            morel_roi = thisMorelAtlas.select_atlas_subset(i);
        
            if num_regions(morel_roi) == 1
                o3 = o2;
                o3.activation_maps = o2.activation_maps(i);
                o3.montage = o2.montage(i);
                morel_roi.montage(o3,'existing_figure','existing_axes', o2.montage{i}.axis_handles,'outline','color',[0,0,0]);
            end
        end
    end
    
    set(gcf,'Tag',[orientation{1},'b']);
end

%% compute overlap statistic

thresh = 0.2;
thisAtlas = myAtlas.resample_space(morel).threshold(0.2);
thisAtlas.probability_maps = [];
thisMorelAtlas = morel;

[newAtlas1, newAtlas2] = deal({});
for i = 1:length(emap)
    newAtlas1{end+1} = thisAtlas.select_atlas_subset(emap{i}{1},'flatten');
    newAtlas2{end+1} = thisMorelAtlas.select_atlas_subset(emap{i}{2},'flatten');
end
thisAtlas = [newAtlas1{:}];
thisMorelAtlas = [newAtlas2{:}];

for i = 1:length(thisAtlas.labels)
    try
        morel_roi = thisMorelAtlas.select_atlas_subset(thisAtlas.labels(i),'exact');
        my_roi = thisAtlas.select_atlas_subset(thisAtlas.labels(i),'exact');
    catch
        continue
    end

    if num_regions(morel_roi) == 1
        coef = 2*sum(morel_roi.dat.*my_roi.dat) / (sum(morel_roi.dat) + sum(my_roi.dat));
        d(i) = coef;
    elseif num_regions(morel_roi) > 1
        error('Too many regions')
    end
end


figure;
hist(d(d~=0))
title('Dice coefficient for LeadDBS Iglesias atlas vs. HCP278 Atlas')
xlabel('Dice Coef.')

