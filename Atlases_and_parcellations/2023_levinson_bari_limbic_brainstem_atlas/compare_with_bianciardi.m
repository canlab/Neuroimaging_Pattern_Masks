clear all;
close all;

addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('/home/bogdan/.matlab/canlab/MasksPrivate'))

limbicAtlas = load_atlas('limbic_brainstem_atlas_fmriprep20').threshold(0.2);

% note, this requires MasksPrivate
biancia = load_atlas('bianciardi_fmriprep20').threshold(0.2);

%% plot overlaps

% equivalence mapping {{iglesias}, {morel}}
emap = {{{'LC_L'},{'L_LC'}},...
    {{'LC_R'},{'R_LC'}},...
    {{'DR'},{'DR_B7'}},...
    {{'PAG'},{'PAG'}},...
    {{'VTA_L'},{'L_VTA_PBP'}},...
    {{'VTA_R'},{'R_VTA_PBP'}},...
    {{'NTS_L'},{'L_VSM'}},...
    {{'NTS_R'},{'R_VSM'}}};

[newAtlas1, newAtlas2] = deal({});
for i = 1:length(emap)
    newAtlas1{end+1} = limbicAtlas.select_atlas_subset(emap{i}{1},'flatten');
    newAtlas2{end+1} = biancia.select_atlas_subset(emap{i}{2},'flatten');
end
limbicBSAtlas = [newAtlas1{:}];
bianciaAtlas = [newAtlas2{:}];

for orientation = {'saggital','coronal','axial'}
    %%
    o2 = limbicBSAtlas.montage('transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(limbicBSAtlas)
        try
            biancia_roi = bianciaAtlas.select_atlas_subset(i).threshold(0.2);
        
            if num_regions(biancia_roi) == 1
                o3 = o2;
                o3.activation_maps = o2.activation_maps(i);
                o3.montage = o2.montage(i);
                biancia_roi.montage(o3,'existing_figure','existing_axes', o2.montage{i}.axis_handles,'outline','color',[0,0,0]);
            end
        end
    end
    
    drawnow()
end
