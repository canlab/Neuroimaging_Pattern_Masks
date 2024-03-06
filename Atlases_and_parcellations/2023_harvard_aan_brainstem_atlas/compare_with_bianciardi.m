clear all;
close all;

addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('/home/bogdan/.matlab/canlab/MasksPrivate'))

aan = load_atlas('harvard_aan_fmriprep20');

% note, this requires MasksPrivate
biancia = load_atlas('bianciardi_fmriprep20');

%% plot overlaps

% equivalence mapping {{iglesias}, {morel}}
emap = {{{'L_LC'},{'L_LC'}},...
    {{'R_LC'},{'R_LC'}},...
    {{'L_LDTg'},{'L_LDTg_CGPn'}},...
    {{'R_LDTg'},{'R_LDTg_CGPn'}},...
    {{'L_PBC'},{'L_MPB','L_LPB'}},...
    {{'R_PBC'},{'R_MPB','R_LPB'}},...
    {{'L_PTg'},{'L_PTg'}},...
    {{'R_PTg'},{'R_PTg'}},...
    {{'L_PnO'},{'L_PnO_PnC_B5'}},...
    {{'R_PnO'},{'R_PnO_PnC_B5'}},...
    {{'L_mRt'},{'L_isRt'}},...
    {{'R_mRt'},{'R_isRt'}},...
    {{'DR'},{'DR_B7'}},...
    {{'MnR'},{'MnR_B6_B8','PMnR_B6_B8'}},...
    {{'PAG'},{'PAG'}},...
    {{'VTA'},{'L_VTA_PBP','R_VTA_PBP'}}};

[newAtlas1, newAtlas2] = deal({});
for i = 1:length(emap)
    newAtlas1{end+1} = aan.select_atlas_subset(emap{i}{1},'flatten');
    newAtlas2{end+1} = biancia.select_atlas_subset(emap{i}{2},'flatten');
end
aanAtlas = [newAtlas1{:}];
bianciaAtlas = [newAtlas2{:}];

for orientation = {'saggital','coronal','axial'}
    %%
    o2 = aanAtlas.montage('transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(aanAtlas)
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
    set(gcf,'Tag',orientation{1});
    
    drawnow()
end
