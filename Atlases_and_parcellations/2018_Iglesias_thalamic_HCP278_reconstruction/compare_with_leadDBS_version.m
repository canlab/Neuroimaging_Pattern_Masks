% J Eugenio Iglesias provided LeadDBS with a copy of the thalamic atlas
% probability space in MNI152NLin2009bAsym space. He described it like so,
% 'I hacked some probabilistic maps for Andy Horn that are in the latest
% version of LeadDBS.'
% I ripped those atlas volumes out of LeadDBS (install it and check in 
% 'INSTALL_ROOT/application/LeadDBS_mcr/LeadDBS/packages/leaddbs/templates/space/MNI152NLin2009bAsym/atlases/Probabilistic
% Thalamus Atlas (Iglesias 2018)') and this script compares those
% segmentations with my own HCP derived ones.

addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))

% path to where I placed the 'Probabilistic Thalamus Atlas (Iglesias 2018)'
% directory
leaddbsRoot='/home/bogdan/MyDocuments/canlab/atlases/_resources/LeadsDBS/';

myAtlas = load_atlas('iglesias_thal_fmriprep20');

ref = fmri_data(which('MNI152NLin2009cAsym_T1_1mm.nii.gz'));

fnames = dir([leaddbsRoot, 'Probablistic Thalamus Atlas (Iglesias 2018)/lh/*nii.gz']);

leadsAtlas = {};
leadsLabel = {};
for i = 1:length(fnames)
    leadsAtlas{end+1} = fmri_data([fnames(i).folder, '/', fnames(i).name]).resample_space(ref);
    leadsLabel{end+1} = ['L_', strrep(strrep(strrep(strrep(fnames(i).name, '.nii.gz',''),'-','_'),'Sg','SG'),'MV(Re)','MV_Re')];
end

fnames = dir([leaddbsRoot, 'Probablistic Thalamus Atlas (Iglesias 2018)/rh/*nii.gz']);

for i = 1:length(fnames)
    leadsAtlas{end+1} = fmri_data([fnames(i).folder, '/', fnames(i).name]).resample_space(ref);
    leadsLabel{end+1} = ['R_', strrep(strrep(strrep(strrep(fnames(i).name, '.nii.gz',''),'-','_'),'Sg','SG'),'MV(Re)','MV_Re')];
end



leadsAtlas = cat(leadsAtlas{:});

% renormalize
p_total = sum(leadsAtlas.dat,2);
leadsAtlas.dat(p_total>1,:) = leadsAtlas.dat(p_total>1,:)./p_total(p_total>1);
leadsAtlas = atlas(leadsAtlas,'labels',leadsLabel);

%% Plot leadsAtlas
cmap = scn_standard_colors(num_regions(leadsAtlas)/2);
cmap = cell2mat(cat(2,cmap'));
cmap = [cmap; cmap];
   
o2 = canlab_results_fmridisplay([], 'full hcp', ...
    'overlay', which(sprintf('%s_T1_1mm.nii.gz',space_description)));
brighten(.6)
    
o2 = montage(atlas2region(leadsAtlas), o2, 'wh_montages', 1:2, 'indexmap', cmap);

%% plot overlap maps

thisAtlas = myAtlas.select_atlas_subset('L_').threshold(0.2);
thisAtlas = thisAtlas.select_atlas_subset(find(~contains(thisAtlas.labels,{'R_L','PuMm','PuMl'})));
thisAtlas.probability_maps = [];
thisLeadsAtlas = leadsAtlas.select_atlas_subset(thisAtlas.labels,'exact').threshold(0.2);
thisLeadsAtlas.probability_maps = [];
for orientation = {'saggital','coronal','axial'}
    %%
    o2 = thisAtlas.montage('transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(thisAtlas)
        try
            leads_roi = thisLeadsAtlas.select_atlas_subset(thisAtlas.labels(i),'exact');
        
            if num_regions(leads_roi) == 1
                o3 = o2;
                o3.activation_maps = o2.activation_maps(i);
                o3.montage = o2.montage(i);
                leads_roi.montage(o3,'existing_figure','existing_axes', o2.montage{i}.axis_handles,'outline','color',[0,0,0]);
            end
        end
    end
    
    drawnow()
end

thisAtlas = myAtlas.select_atlas_subset('R_').threshold(0.2);
thisAtlas = thisAtlas.select_atlas_subset(find(~contains(thisAtlas.labels,{'PuMm','PuMl'})));
thisAtlas.probability_maps = [];
thisLeadsAtlas = leadsAtlas.select_atlas_subset(thisAtlas.labels,'exact').threshold(0.2);
thisLeadsAtlas.probability_maps = [];
for orientation = {'saggital','coronal','axial'}
    %%
    o2 = thisAtlas.montage('transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(thisAtlas)
        try
            leads_roi = thisLeadsAtlas.select_atlas_subset(thisAtlas.labels(i),'exact');
        
            if num_regions(leads_roi) == 1
                o3 = o2;
                o3.activation_maps = o2.activation_maps(i);
                o3.montage = o2.montage(i);
                leads_roi.montage(o3,'existing_figure','existing_axes', o2.montage{i}.axis_handles,'outline','color',[0,0,0]);
            end
        end
    end
    drawnow()
end


%% overlap statistics
[d,r] = deal(zeros(num_regions(myAtlas),1));
for i = 1:length(myAtlas.labels)
    try
    leads_roi = leadsAtlas.select_atlas_subset(myAtlas.labels(i),'exact');
    my_roi = myAtlas.select_atlas_subset(myAtlas.labels(i),'exact');
    catch
        continue
    end

    if num_regions(leads_roi) == 1
        joint_vx = logical(leads_roi.dat | my_roi.dat);
        r(i) = corr(leads_roi.probability_maps(joint_vx), my_roi.probability_maps(joint_vx)); 
    elseif num_regions(leads_roi) > 1
        error('Too many regions')
    end
end


thresh = 0.2;
leadsAtlas = leadsAtlas.threshold(thresh);
leadsAtlas.probability_maps = [];
myAtlas = myAtlas.threshold(thresh);
myAtlas.probability_maps = [];
for i = 1:length(myAtlas.labels)
    try
    leads_roi = leadsAtlas.select_atlas_subset(myAtlas.labels(i),'exact');
    my_roi = myAtlas.select_atlas_subset(myAtlas.labels(i),'exact');
    catch
        continue
    end

    if num_regions(leads_roi) == 1
        coef = 2*sum(leads_roi.dat.*my_roi.dat) / (sum(leads_roi.dat) + sum(my_roi.dat));
        d(i) = coef;
    elseif num_regions(leads_roi) > 1
        error('Too many regions')
    end
end


figure;
hist(d(d~=0))
title('Dice coefficient for LeadDBS Iglesias atlas vs. HCP278 Atlas')
xlabel('Dice Coef.')

figure;
hist(r(r~=0))
title({'Correnation coefficient for LeadDBS Iglesias atlas','vs. HCP278 Atlas'' probabilities'});
xlabel('Pearson r')