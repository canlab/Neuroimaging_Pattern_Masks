% This script was adapted from the single_trial_database prep1 script.
% It requires the CANlab single trial database, with specific setup code.
% Run after running a_set_up_paths... and b_load_check_prep_for_analysis.
% This was designed to extract data for characterizing Tor Wager's 2019
% pain pathways identification specifically.

% Goals:
% - Create and save a table of all subject-level data, including .mat files
% with single-trial data
% -
% Extract single-trial values for each region in pain pathways objects
% Define 1 variable storing data:
%
%   ST_data_cPDM
% --------------------------------------------------------

load pain_pathways_region_obj_with_local_patterns pain_regions_cpdm % has local weights extracted

% Very important to load these things, as we'll use gwcsf add to ST_cleaned here:

cd /Users/torwager/Documents/GitHub/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2019_Wager_pain_pathways/data
load data/pain_pathways_single_trial_data_2019 ST_data_gray_white_csf ST_cleaned


%% SET UP ALL STUDIES

% EDIT HERE TO ADD NEW STUDY
% ------------------------------------------------------------------------
studynames = {'NSF' 'BMRK3' 'BMRK4' 'IE' 'SCEBL' 'ILCP' 'bmrk5_painsound' 'EXP' 'RCP' 'IE2' 'placebo_value_stephan' 'REMI' 'levoderm' };  % LEVO, Romanticpain, BMRK5

studydirs = {'nsf' ...
    'Tor_bmrk3_datashare' ...
    'bmrk4_savedbywani' ...   'bmrk4_pain_ST_SmoothedWithBasis_handfoot' ... % alternative? noted by tor, 7/2019
    'ie_for_tor' ...
    'SCEBL_single_trial_Leonie' ...
    'ILCP_wani' ...
    'bmrk5_painsound' ...
    'Expectancy' ...
    'romantic_pain' ...
    'IE2_NEW' ...
    'placebo_value_stephan' ...
    'REMI_wani' ...
    'levoderm' ...
    };


studymetadatanames = {'meta.mat' ...
    'bmrk3_single_trial_model.mat' ...
    'bmrk4_single_trial_model.mat' ...
    'ie_model_with_vif.mat' ...
    'SCEBLdata_forTor_N26.mat' ...
    'ilcp_metainfo_waniupdate.mat' ...
    'single_trial_bmrk5_meta' ...
    'exp_meta.mat' ...
    'romantic_pain.mat' ...
    'ie2_metadata' ...
    'placebo_value_meta_tor.mat' ...
    'remi_st_metadata.mat' ...
    'levoderm.mat' ...
    };
% ------------------------------------------------------------------------
% NOT INCLUDED
% 'bmrk5_mod' ... 'bmrk5_modulationtasks' ...'bmrk5_mod_meta.mat' ...
% 'APP-fMRI'  'APP-fMRI' 'mediation_variables2.mat'

nstudies = length(studynames);

%% CREATE CANLAB_DATASET OBJECT FOR EACH STUDY
% see plugin script for more info

clear matfiles wh_study

for i = 1:nstudies
    
    studydatapath{i} = [fullfile(datadir, studydirs{i}) filesep];
    
    tmp  = load(fullfile(studydatapath{i}, studymetadatanames{i}));
    n = fieldnames(tmp);
    studymetadata{i} = tmp.(n{1});
    
    studysubjects{i} = studymetadata{i}.subjects;
    
    wh_study{i, 1} = i * ones(length(studysubjects{i}), 1);  % study index number for each subject
    
    tmp = strvcat(studymetadata{i}.dat_obj{:});
    matfiles{i} = [repmat(studydatapath{i}, size(tmp, 1), 1) tmp ];
    
    studymetadata{i}.studyname = studynames(i * ones(size(wh_study{i}, 1), 1));
    
    % AD-HOC FIXES
    if ~isfield(studymetadata{i}, 'images') && isfield(studymetadata{i}, 'imgs')
        studymetadata{i}.images = studymetadata{i}.imgs;
    end
    
    if ~isfield(studymetadata{i}, 'vifs') && isfield(studymetadata{i}, 'vif')
        studymetadata{i}.vifs = studymetadata{i}.vif;
    end
    
    % ------------------------------------
    meta_to_canlab_dataset_plugin;
    % ------------------------------------
    
    study_canlab_dataset{i} = DAT;
    
    %save(fullfile(studydatapath{i}, 'canlab_dataset_object.mat'), 'DAT');
end

% This script sets up the canlab_dataset object for the data in this folder.

%save(fullfile(resultsdir, 'canlab_datasets_and_metadata.mat'), 'studymetadata', 'study_canlab_dataset')

%% Get file names for each subject's mat file to load

dorestart = true;

if dorestart
    % Initialize variables
    all_data_matfiles = {};
    [all_study_indx, all_subj_indx, all_data_matfile_ok] = deal([]);
    
    [ST_data_cPDM] = deal({});
    
    indx = 1;
end

for i = 1:nstudies
    
    matfile = get_var(study_canlab_dataset{i}, 'ST_matfile');
    
    for j = 1:length(matfile)
        
        fprintf('Study %3.0f %s Subj %3.0f\n', i, studynames{i}, j);
        
        % save data for table
        % ---------------------------------------
        all_data_matfiles{end+1} = matfile{j};
        all_study_indx(end+1) = i;
        all_subj_indx(end+1) = indx;
        
        all_data_matfile_ok(end+1) = false;
        % ---------------------------------------
        
        clear dat
        
        if exist(matfile{j}, 'file')
            
            all_data_matfile_ok(end) = true;
            
            dat = load(matfile{j});
            n = fieldnames(dat);
            dat = dat.(n{1});
            
            % Extract values for each region in painpathways
            %
            % Define  variables:
            %   ST_data_cPDM
            % --------------------------------------------------------
            tic
            [regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_cpdm, dat);
            
            ST_data_region_averages{indx} = cat(2, regions_with_testdata.dat);  % trials x regions
            ST_data_cPDM{indx} = cat(2, local_pattern_responses{:});  % trials x regions
            
            toc
            
            % kludge - stephan value study #38
            if i == 11 && j == 38
                tmp = study_canlab_dataset{i}.Event_Level.data{j};
                [wasnan] = nanremove(tmp(:, 1));
                tmp(wasnan, :) = [];
                
                %                 ST_data_region_averages{indx}(wasnan, :) = [];
                %                 ST_data_PDM1_response{indx}(wasnan, :) = [];
                %                 ST_data_NPS_response{indx}(wasnan, :) = [];
                %                 ST_data_SIIPS_response{indx}(wasnan, :) = [];
                %                 ST_data_painpathways_finegrained{indx}(wasnan, :) = [];
                
                study_canlab_dataset{i}.Event_Level.data{j} = tmp;
                
            end
            
            
        else
            disp('matfile missing');
            
            ST_data_cPDM{indx} = [];
            
            
        end
        
        toc
        
        indx = indx + 1;
        
    end
    
end %

%% Make table

Matfile = all_data_matfiles';
Study = all_study_indx';
Subject = all_subj_indx';
Mat_OK = all_data_matfile_ok';

Subject_Table = table(Study, Subject, Mat_OK, Matfile);

Studyname = studynames(Subject_Table.Study)';
Subject_Table.Studyname = Studyname;

%% Save

save pain_pathways_single_trial_data_2019 -append ST_data_cPDM

%% These are on different scales.  Rescale by z-scoring each variable
% and concatenate.

% Gray/white/CSF - to regress out WM and CSF averages. Z-score within
% subject.
mydata = cellfun(@zscore, ST_data_gray_white_csf, 'UniformOutput', false);
all_ST_data_gray_white_csf_Z = cat(1, mydata{:});
%figure; imagesc(zscore(all_ST_data_gray_white_csf_Z))
title('Gray - white - CSF averages')

% Regress out white-matter and CSF averages (doens't make much difference)
X = all_ST_data_gray_white_csf_Z(:, 2:3);
resid_fcn = @(Y) Y - X * pinv(X) * Y;


%% cPDM
% -----------------------------------------------------
% Z-score within-subject

mydata = cellfun(@windsorize_matrix_columnwise, ST_data_cPDM, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cellfun(@zscore, mydata, 'UniformOutput', false); % Z-score within subject
mydata = cat(1, mydata{:});

% Regress out white-matter and CSF averages (doens't make much difference)
v = var(mydata);
mydata = resid_fcn(mydata);
vr = var(mydata);
vexplained = 1 - vr ./ v;

% Windsorize again to Z = 3 (within-subject, because Z-scores are calculated w/i)
wh_out = mydata > 3; percent_above = 100 * sum(wh_out(:)) ./ prod(size(mydata));
mydata(wh_out) = 3;
wh_out = mydata < -3; percent_below = 100 * sum(wh_out(:)) ./ prod(size(mydata));
mydata(wh_out) = -3;
descrip = sprintf('windsorized %3.0f%% of trials to 3 standard deviations\n', percent_above + percent_below);
disp(descrip)

ST_cleaned.cPDM = mydata;
ST_cleaned.cPDM_history = sprintf('Z-score trials within subject column-wise\n%s\nRegressed out WM+CSF\nexplained %3.0f%% of variance on average', descrip, 100*mean(vexplained) );
ST_cleaned.cPDM_history = textwrap({ST_cleaned.cPDM_history}, 80);

figure; imagesc(mydata); colorbar
title('cPDM')

% -----------------------------------------------------


%% Add labels

ST_cleaned.cpdm_labels = {pain_regions_cpdm.shorttitle};

%% Remove bad trials (sound BMRK5)

ST_cleaned.cPDM(ST_cleaned.wasremoved, :) = [];
%% Save

save data/pain_pathways_single_trial_data_2019 -append ST_cleaned





