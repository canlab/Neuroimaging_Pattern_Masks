function [meta, regional_data] = preparation_extract_data(study_canlab_dataset, studymetadata)
% get the meta information including temp, OK trials, GM mean, white mean and CSF values
% get regional average activity

%% Extract meta data
%--------------------------------------------------------------------------

Nstudies = length(study_canlab_dataset);
extract_idx = [1 2 22 23 24 43]; % index of features need to extract 
extract_name = {'ratings' 'temp' 'GM mean' 'WM mean' 'CSF mean' 'OK trials'}; % name of features need to extract 

study_name = {'NSF' 'bmrk3' 'bmrk4' 'IE' 'SCEBL' 'ILCP' 'bmrk5' 'EXP'...
    'RCP' 'IE2' 'placebo_stephan' 'REMI' 'levoderm'};

rm_idx = [9, 11, 13]; % remove RCP, placebo_stephan, levoderm
% rm_idx = [3, 4, 5, 7, 9, 10, 11, 13]; % index of removed studies
study_name(rm_idx) = [];
study_name(cellfun(@isempty, study_name)) = [];

for i = 1:Nstudies %find which study has placebo/drug trials
    
    Nsubjets = length(study_canlab_dataset{i}.Event_Level.data);
    Subj_id{i} = studymetadata{i}.dat_obj;
    
    if isfield(studymetadata{i}, 'placebotrial')
        
        placebo_drug{i} = studymetadata{i}.placebotrial;
        
    end
       
    if isfield(studymetadata{i}, 'drug')
        
        placebo_drug{i} = studymetadata{i}.drug;
        
    elseif isfield(studymetadata{i}, 'placebo')
        
        placebo_drug{i} = studymetadata{i}.placebo;
        
    else
        
        placebo_drug{i} = [];
    
    end
    
    for j = 1:Nsubjets
        
        data_tmp = study_canlab_dataset{i}.Event_Level.data{j};
        extract_data_tmp = data_tmp(:, extract_idx);
        meta_data{i}{j} = extract_data_tmp;
        
    end
     
end

Subj_id(rm_idx) = []; Subj_id(cellfun(@isempty, Subj_id)) = [];
placebo_drug(rm_idx) = []; placebo_drug(cellfun(@isempty, placebo_drug)) = [];
meta_data(rm_idx) = []; meta_data(cellfun(@isempty, meta_data)) = [];

meta.study_name = study_name;
meta.var_name = extract_name;
meta.Subj_id = Subj_id;
meta.meta_data = meta_data;
meta.placebo_drug_trial{4} = placebo_drug{1};
meta.placebo_drug_trial{10} = placebo_drug{2};

%% Extract regional mean data
%--------------------------------------------------------------------------

dir = 'E:\Matlab_stuff\workspace\data\pain\datasets';
template = 'E:\Matlab_stuff\workspace\template\BN_atlas\BN_Atlas_274\BN_Atlas_274_noCb_uint16.nii';

for i = 1:length(study_name)
    
    for j = 1:length(Subj_id{i})
        
        subj_dir = strcat([dir '\' study_name{i} '\' Subj_id{i}{j}]);
        load(subj_dir);
        
        if strcmp('SCEBL', study_name{i}) == 1
            
           [~, roi_val] = canlab_connectivity_preproc(painbetas, 'extract_roi', template, 'no_preproc');
           
        elseif strcmp('REMI', study_name{i}) == 1
            
           [~, roi_val] = canlab_connectivity_preproc(dat_obj, 'extract_roi', template, 'no_preproc');
            
        else
            
           [~, roi_val] = canlab_connectivity_preproc(dat, 'extract_roi', template, 'no_preproc');
            
        end
        
        roi_tmp = roi_val{1}.dat; 
        
        regional_data{i}{j} = roi_tmp;
        clear roi_tmp;
        
    end

end







