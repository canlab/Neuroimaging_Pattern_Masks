% create_brainstem_atlas_group
% 
% Create a brain atlas with anatomically defined region groups, from
% various atlases/papers.  Uses canlab_load_ROI
%
% Notes: functional atlases probably do not [yet] have very good
% subdivisions, and there is a clear demarcation of functions, inputs, and
% outputs by anatomical subnuclei.

% Define: 1 mm space by default, based on HCP image
% This initial image covers the whole space

bstem_atlas = atlas(which('canlab_brainstem.img'));
bstem_atlas = threshold(bstem_atlas, .2);

bstemimg = fmri_data(which('brainstem_mask_tight_2018.img'));
bstem_atlas = apply_mask(bstem_atlas, bstemimg);

% this has some other regions in it (e.g., mammillary bodies), so would be better to use the
% 'noreplace' option when merging it with other atlases.

orthviews(bstem_atlas);

%see also: bstemimg.fullpath = fullfile(pwd, 'brainstem_mask_tight_2018.img');

%% Shen regions: Fill in parcels not assigned to a named region

shen = load_atlas('shen');

shen = apply_mask(shen, bstemimg);

r = atlas2region(shen);
orthviews(shen)
%[r, labels] = cluster_names(r,true);

labels1 = {'xxx' 'Shen_Midb_Rrv' 'Shen_CerPed_R' 'Shen_CerPed_R' 'Shen_CerPed_R' 'xxx' 'xxx' 'xxx'};  % 1:8
%for i = 1:8, orthviews(r(i), {[1 0 0]}); labels1{i}, input(''), end

labels2 = {'xxx' 'Shen_Midb_Rd' 'xxx' 'Shen_Med_R' 'Shen_Pons_R' 'Shen_Pons_Rcv' 'Shen_Midb_R' 'Shen_Pons_Rcd'};
rr = r(9:16);
%for i = 1:8, orthviews(rr(i), {[1 0 0]}); labels2{i}, input(''), end

labels3 = {'xxx' 'xxx' 'Shen_CerPed_L' 'xxx' 'xxx' 'xxx' 'xxx' 'Shen_Pons_Lrd'};
rr = r(17:24);
%for i = 1:8, orthviews(rr(i), {[1 0 0]}); labels3{i}, input(''), end

labels4 = {'xxx' 'xxx' 'xxx' 'Shen_Midb_Ld' 'xxx' 'Shen_Midb_L' 'Shen_Med_L' 'Shen_Pons_Lcd' 'Shen_Pons_L'};
rr = r(25:33);
%for i = 1:9, orthviews(rr(i), {[1 0 0]}); labels4{i}, input(''), end

labels = [labels1 labels2 labels3 labels4];
wh_omit = strcmp(labels, 'xxx');

r(wh_omit) = [];
labels(wh_omit) = [];

shen = remove_atlas_region(shen, find(wh_omit));
shen.labels_3 = shen.labels_2;
shen.labels_2 = shen.labels;
shen.labels = labels;

% to find clusters by hand:
% [~,wh] = find_closest_cluster(r, spm_orthviews('pos')) 

% Reorder Shen regions, leaving out peduncles

[~, left] = select_atlas_subset(shen, {'_L'});
[~, right] = select_atlas_subset(shen, {'_R'});

[~, wh_midb] = select_atlas_subset(shen, {'Midb'});
wh_midb = [find([wh_midb & left]) find([wh_midb & right])];

[~, wh_pons] = select_atlas_subset(shen, {'Pons'});
wh_pons = [find([wh_pons & left]) find([wh_pons & right])];

[~, wh_med] = select_atlas_subset(shen, {'Med'});
wh_med = [find([wh_med & left]) find([wh_med & right])];


wh_order = [wh_midb wh_pons wh_med];

shen = reorder_atlas_regions(shen, wh_order);


%% Add shen brainstem, replacing old ones

bstem_atlas = merge_atlases(bstem_atlas, shen, 'always_replace');


%% add other regions
% ...replacing voxels where new one overlaps

% also include regions in other atlases that we want to remove here - so
% that we remove these voxels

% to-do: 'pbn' 'nts'

regionnames = {'pag' 'sc' 'ic' 'drn' 'mrn' 'PBP' 'sn' 'VTA' 'rn'  'pbn' 'lc' 'rvm' 'nrm' 'dmnx_nts' 'ncf' 'ncs_B6_B8' 'nrp_B5' 'nuc_ambiguus' 'medullary_raphe' 'spinal_trigeminal'};

% NEW ONES TOO

for i = 1:length(regionnames)
    regionname = regionnames{i};
    
    [~, roi_atlas] = canlab_load_ROI(regionname);
    orthviews(roi_atlas);

    bstem_atlas = merge_atlases(bstem_atlas, roi_atlas, 'always_replace');
    
end

% Fix - not sure why some labels not saving
% bstem_atlas.labels(end-2:end) = {'R_LC' 'L_LC' 'rvm'};

atlas_obj = bstem_atlas;

%% Adjust labels
% make more consistent with other atlases
% relabel L and R

pat = 'Reg_1';
atlas_obj.labels = regexprep(atlas_obj.labels, pat, 'other');

pat = 'Shen_';
atlas_obj.labels = regexprep(atlas_obj.labels, pat, '');

atlas_obj = atlas_add_L_R_to_labels(atlas_obj);


%% Add references

references = {'Shen, X., F. Tokoglu, X. Papademetris, and R. T. Constable. 2013. ?Groupwise Whole-Brain Parcellation from Resting-State fMRI Data for Network Node Identification.? NeuroImage 82 (November): 403?15.'
    'Pauli, Wolfgang M., Amanda N. Nili, and J. Michael Tyszka. 2018. ?A High-Resolution Probabilistic in Vivo Atlas of Human Subcortical Brain Nuclei.? Scientific Data 5 (April): 180063.'
    'Fairhurst, Merle, Katja Wiech, Paul Dunckley, and Irene Tracey. 2007. ?Anticipatory Brainstem Activity Predicts Neural Processing of Pain in Humans.? Pain 128 (1-2):101?10.'
    'Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.'
    'Zambreanu, L., R. G. Wise, J. C. W. Brooks, G. D. Iannetti, and I. Tracey. 2005. ?A Role for the Brainstem in Central Sensitisation in Humans. Evidence from Functional Magnetic Resonance Imaging.? Pain 114 (3):397?407.'
    'Keuken, M. C., P-L Bazin, L. Crown, J. Hootsmans, A. Laufer, C. Müller-Axt, R. Sier, et al. 2014. ?Quantifying Inter-Individual Anatomical Variability in the Subcortex Using 7 T Structural MRI.? NeuroImage 94 (July): 40?46.'
    'Beliveau, Vincent, Claus Svarer, Vibe G. Frokjaer, Gitte M. Knudsen, Douglas N. Greve, and Patrick M. Fisher. 2015. ?Functional Connectivity of the Dorsal and Median Raphe Nuclei at Rest.? NeuroImage 116 (August): 187?95.'
    'Sclocco, Roberta, Florian Beissner, Gaelle Desbordes, Jonathan R. Polimeni, Lawrence L. Wald, Norman W. Kettner, Jieun Kim, et al. 2016. ?Neuroimaging Brainstem Circuitry Supporting Cardiovagal Response to Pain: A Combined Heart Rate Variability/ultrahigh-Field (7 T) Functional Magnetic Resonance Imaging Study.? Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences 374 (2067). rsta.royalsocietypublishing.org. https://doi.org/10.1098/rsta.2015.0189.'
    'Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.'
    'Keren, Noam I., Carl T. Lozar, Kelly C. Harris, Paul S. Morgan, and Mark A. Eckert. 2009. ?In Vivo Mapping of the Human Locus Coeruleus.? NeuroImage 47 (4): 1261?67.'
    };

atlas_obj.references = char(references);

      
% 'pag'        Periaqueductal gray, hand-drawn (Tor Wager 2018, mask out aqueduct/Keuken2014)
% 'sc'         Superior colliculus, hand-drawn (Tor Wager 2018, mask out aqueduct/Keuken2014)
% 'ic'         Inferior colliculus, hand-drawn (Tor Wager 2018, mask out aqueduct/Keuken2014)
% 'drn'        Dorsal raphe nucleus, coords from Beliveau, 2015. mask out aqueduct/Keuken2014 ?Functional Connectivity of the Dorsal and Median Raphe Nuclei at Rest.? NeuroImage 116 (August). Elsevier:187?95.
% 'mrn'        Median raphe nucleus, coords from Beliveau, 2015. mask out aqueduct/Keuken2014 ?Functional Connectivity of the Dorsal and Median Raphe Nuclei at Rest.? NeuroImage 116 (August). Elsevier:187?95.
% 'PBP'        Parabrachial pigmented nuc.      % Pauli 2017 BioArxiv subcortical atlas
% 'sn'         Substantia Nigra; Keuken 2014
% 'SNc'        Substantia Nigra compacta        % Pauli 2017 BioArxiv subcortical atlas
% 'SNr'        Substantia Nigra reticularis     % Pauli 2017 BioArxiv subcortical atlas
% 'VTA'        Ventral tegmental area           % Pauli 2017 BioArxiv subcortical atlas
% 'rn'         Red nucleus; Keuken 2014
% 'pbn'        Parabrachial complex; Fairhurst, Merle, Katja Wiech, Paul Dunckley, and Irene Tracey. 2007. ?Anticipatory Brainstem Activity Predicts Neural Processing of Pain in Humans.? Pain 128 (1-2):101?10.
% 'lc'         Locus coeruleus; Keren 2009, 2SD image
% 'rvm_old'    Hand-drawn rostral ventral medulla (Tor) in anatomical rvm
% 'rvm'        Rostral ventral medulla from Brooks et al. 2016(??)
% 'nts'        Nuc. tractus solitarius (rough; hand-drawn, Tor)
% 'olive'      Inferior olive; MISSING
% 'nrm'        Nuc. raphe magnus; % Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% 'ncf'        Nuc. cuneiformis; Zambreanu, L., R. G. Wise, J. C. W. Brooks, G. D. Iannetti, and I. Tracey. 2005. ?A Role for the Brainstem in Central Sensitisation in Humans. Evidence from Functional Magnetic Resonance Imaging.? Pain 114 (3):397?407.
% 'ncs_B6_B8'  Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% 'nrp_B5'     Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% 'nuc_ambiguus' Sclocco, Roberta, Florian Beissner, Gaelle Desbordes, Jonathan R. Polimeni, Lawrence L. Wald, Norman W. Kettner, Jieun Kim, et al. 2016. ?Neuroimaging Brainstem Circuitry Supporting Cardiovagal Response to Pain: A Combined Heart Rate Variability/ultrahigh-Field (7 T) Functional Magnetic Resonance Imaging Study.? Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences 374 (2067). rsta.royalsocietypublishing.org. https://doi.org/10.1098/rsta.2015.0189.
% 'dmnx_nts'    Sclocco, Roberta, Florian Beissner, Gaelle Desbordes, Jonathan R. Polimeni, Lawrence L. Wald, Norman W. Kettner, Jieun Kim, et al. 2016. ?Neuroimaging Brainstem Circuitry Supporting Cardiovagal Response to Pain: A Combined Heart Rate Variability/ultrahigh-Field (7 T) Functional Magnetic Resonance Imaging Study.? Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences 374 (2067). rsta.royalsocietypublishing.org. https://doi.org/10.1098/rsta.2015.0189.
% 'vep'          % Pauli 2017 BioArxiv CIT168 subcortical atlas 
% 'medullary_raphe' Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.
% 'spinal_trigeminal' Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.
%

%% REMOVE general region - not needed

atlas_obj = remove_atlas_region(atlas_obj, {'other'});

%% Save dir

savedir = what('2018_Wager_combined_atlas');
savedir = savedir.path;

cd(savedir)

%% save object

atlas_name = 'brainstem_combined';

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save(savename, 'atlas_obj');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

r = atlas2region(atlas_obj);
labels = atlas_obj.labels;

for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    set(han,'FaceAlpha', .5)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end

 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'multirow', 1);
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

end
 