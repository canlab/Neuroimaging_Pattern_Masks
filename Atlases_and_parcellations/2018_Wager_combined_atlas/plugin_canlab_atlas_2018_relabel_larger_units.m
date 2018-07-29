%% PLUGIN: LABEL CANLAB2018 ATLAS LABELS
% ADD larger organziational labels and include them as labels_2.
% These can be used in tables to report additional fields or sort by them.
%
% Load, relabel,and save

% tor wager, July 2018
% this is a one-time script used in atlas creation, in
% script_2018_Wager_combined_atlas.  

%% Save dir

savedir = what('2018_Wager_combined_atlas');
savedir = savedir.path;

cd(savedir)

%% LOAD combined atlas

load(fullfile(savedir, 'CANlab_combined_atlas_object_2018.mat'), 'atlas_obj');

new_labels = relabel_larger_units(atlas_obj);

atlas_obj.labels_2 = new_labels;

% Replace "Cortex" with labels from Yeo 17-network (16 unique) atlas:
atlas_obj = relabel_cortex_with_functional_networks(atlas_obj);

% replace label_descriptions so these will show up in tables, in region.table
atlas_obj.labels_3 = atlas_obj.label_descriptions';
atlas_obj.label_descriptions = atlas_obj.labels_2';

%% SAVE combined atlas

save(fullfile(savedir, 'CANlab_combined_atlas_object_2018.mat'), 'atlas_obj');

atlas_obj.fullpath = fullfile(pwd, 'CANlab_2018_combined_atlas.nii');
write(atlas_obj);
%% LOAD combined atlas - downsampled (used in tables)

load(fullfile(savedir, 'CANlab_combined_atlas_object_2018_2mm.mat'), 'atlas_obj');

new_labels = relabel_larger_units(atlas_obj);

atlas_obj.labels_2 = new_labels;

% Replace "Cortex" with labels from Yeo 17-network (16 unique) atlas:
atlas_obj = relabel_cortex_with_functional_networks(atlas_obj);

% replace label_descriptions so these will show up in tables, in region.table
atlas_obj.labels_3 = atlas_obj.label_descriptions';
atlas_obj.label_descriptions = atlas_obj.labels_2';

%% SAVE combined atlas

save(fullfile(savedir, 'CANlab_combined_atlas_object_2018_2mm.mat'), 'atlas_obj');

atlas_obj.fullpath = fullfile(pwd, 'CANlab_2018_combined_atlas_2mm.nii');
write(atlas_obj);

%%
function new_labels = relabel_larger_units(atlas_obj)

new_labels = cell(size(atlas_obj.labels));

% Labels_2 has broader units, for sorting table rows

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Ctx')); 
new_labels(wh) = {'Cortex'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Cblm'));
new_labels(wh) = {'Cerebellum'};

% BG
% ----------------------------------------------------------------------
wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Caudate'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Putamen'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'V_Striatum'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'NAC_'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'GPe'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'GPi'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'STN'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'VeP'));
new_labels(wh) = {'Basal_ganglia'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Cau_')); % redundant? should be removed or merged in future probably
new_labels(wh) = {'Basal_ganglia'};

% Other telencephalon
% ----------------------------------------------------------------------
wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Amygdala'));
new_labels(wh) = {'Amygdala'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Hippocampus'));
new_labels(wh) = {'Hippocampus'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Subic'));
new_labels(wh) = {'Hippocampus'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'BST_'));
new_labels(wh) = {'Basal_forebrain'};



% Thalamus and other diencephalon
% ----------------------------------------------------------------------
wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Thal'));
new_labels(wh) = {'Diencephalon'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Haben'));
new_labels(wh) = {'Diencephalon'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Hythal')); % should already be labeled but including again for clarity
new_labels(wh) = {'Diencephalon'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Mamm_')); % should already be labeled but including again for clarity
new_labels(wh) = {'Diencephalon'};



% Brainstem
% ----------------------------------------------------------------------

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Bstem'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'raphe'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'PBP_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'sn_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'pbn_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'lc_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'rvm_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'dmnx_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'ambiguus'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'spinal_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'IC_')); % inf collic - should have been labeled Bstem before but missed
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'ncf_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'ncs_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'nrp_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'nrm_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'VTA_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'rn_'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'nrm'));
new_labels(wh) = {'Brainstem'};

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Ctx')); % do this again because IC_ screws up one of the cortical labels
new_labels(wh) = {'Cortex'};

end % function




function atlas_obj = relabel_cortex_with_functional_networks(atlas_obj)


%% RELABEL cortical regions by network

% Load atlas with network labels: Use Yeo 17 network labels (16 unique)

atlas_to_add = load_atlas('yeo17networks');

% Match each original region in atlas_obj to the best-matching new label:

r = atlas2region(atlas_obj);
[r, region_table, table_legend_text] = autolabel_regions_using_atlas(r, atlas_to_add);

new_labels_raw = {r.shorttitle};  % network labels with original label text

% Re-label cortical regions only 

wh = ~cellfun(@isempty, strfind(atlas_obj.labels, 'Ctx')); 

%%
new_labels = new_labels_raw(wh);
new_labels = strrep(new_labels, 'LH_', 'Cortex_'); % remove hemisphere labels - they are already in individual regions
new_labels = strrep(new_labels, 'RH_', 'Cortex_'); % add Cortex for context.

% replace with more descriptive names - also focus on anatomical rather
% than functional labels, e.g., "Fronto-parietal" instead of "control"
% though not always possible with current labels in literature...

new_labels = strrep(new_labels, 'Default', 'Default_Mode');
new_labels = strrep(new_labels, 'Cont', 'Fronto_Parietal');
new_labels = strrep(new_labels, 'SalVentAttn', 'Ventral_Attention');
new_labels = strrep(new_labels, 'DorsAttn', 'Dorsal_Attention');
new_labels = strrep(new_labels, 'SomMot', 'Somatomotor');
new_labels = strrep(new_labels, 'VisCent', 'Visual_Central');
new_labels = strrep(new_labels, 'VisPeri', 'Visual_Peripheral');
new_labels = strrep(new_labels, 'TempPar', 'Temporal_Parietal');

atlas_obj.labels_2(wh) = new_labels;


end % subfunction


