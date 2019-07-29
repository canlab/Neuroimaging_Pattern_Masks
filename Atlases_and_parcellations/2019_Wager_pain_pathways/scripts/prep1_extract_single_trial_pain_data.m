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
% Define 5 variables stoing data:
%
%   ST_data_region_averages
%   ST_data_PDM1_response
%   ST_data_NPS_response
%   ST_data_SIIPS_response
%   ST_data_painpathways_finegrained
%
% ST_data_pain_temp contains pain ratings (raw, non-normalized) and temperatures for each trial
%
% ST_data_gray_white_csf                cell array (1 per subject) with ST [gray white CSF] means
% ST_data_gray_white_csf_components     cell array (1 per subject) with ST[gray x 5 white x 5 CSF x 5] principal components
% --------------------------------------------------------


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

dorestart = false;

if dorestart
    % Initialize variables
    all_data_matfiles = {};
    [all_study_indx, all_subj_indx, all_data_matfile_ok] = deal([]);
    
    [ST_data_region_averages, ST_data_PDM1_response, ST_data_NPS_response, ST_data_SIIPS_response, ST_data_painpathways_finegrained, ST_data_pain_temp, ST_data_gray_white_csf, ST_data_gray_white_csf_components] = deal({});
    
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
            % Define 5 variables:
            %   ST_data_region_averages
            %   ST_data_PDM1_response
            %   ST_data_NPS_response
            %   ST_data_SIIPS_response
            %   ST_data_painpathways_finegrained
            % --------------------------------------------------------
            tic
            [regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_pdm1, dat);
            
            ST_data_region_averages{indx} = cat(2, regions_with_testdata.dat);  % trials x regions
            ST_data_PDM1_response{indx} = cat(2, local_pattern_responses{:});  % trials x regions
            
            [~, local_pattern_responses] = extract_data(pain_regions_nps, dat);
            ST_data_NPS_response{indx} = cat(2, local_pattern_responses{:});  % trials x regions
            
            [~, local_pattern_responses] = extract_data(pain_regions_siips, dat);
            ST_data_SIIPS_response{indx} = cat(2, local_pattern_responses{:});  % trials x regions
            
            regions_with_testdata = extract_roi_averages(dat, pain_pathways_finegrained);
            ST_data_painpathways_finegrained{indx} = cat(2, regions_with_testdata.dat);  % trials x regions
            
            tmp = study_canlab_dataset{i}.Event_Level.data{j}; % trial-level pain and temperature
            ST_data_pain_temp{indx} = tmp(:, 1:2);
            
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
            
            % Extract gray, white, CSF components
            % ---------------------------------------------------------------
            [gray_white_csf, gwc_components] = extract_gray_white_csf(dat);
            
            ST_data_gray_white_csf{indx} = gray_white_csf;
            ST_data_gray_white_csf_components{indx} = cat(2, gwc_components{:});
            
        else
            disp('matfile missing');
            
            ST_data_region_averages{indx} = [];
            ST_data_PDM1_response{indx} = [];
            ST_data_NPS_response{indx} = [];
            ST_data_SIIPS_response{indx} = [];
            ST_data_painpathways_finegrained{indx} = [];
            
            ST_data_gray_white_csf{indx} = [];
            ST_data_gray_white_csf_components{indx} = [];
            
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

cd /Users/torwager/Documents/GitHub/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2019_Wager_pain_pathways/data

save pain_pathways_single_trial_data_2019 Subject_Table ST_*

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


%% NPS
% -----------------------------------------------------
% Z-score within-subject

mydata = cellfun(@windsorize_matrix_columnwise, ST_data_NPS_response, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
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

ST_cleaned.NPS = mydata;
ST_cleaned.NPS_history = sprintf('Z-score trials within subject column-wise\n%s\nRegressed out WM+CSF\nexplained %3.0f%% of variance on average', descrip, 100*mean(vexplained) );
ST_cleaned.NPS_history = textwrap({ST_cleaned.NPS_history}, 80);

figure; imagesc(mydata); colorbar
title('NPS')

% -----------------------------------------------------

%% SIIPS
% -----------------------------------------------------
% Z-score within-subject

mydata = cellfun(@windsorize_matrix_columnwise, ST_data_SIIPS_response, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
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

ST_cleaned.SIIPS = mydata;
ST_cleaned.SIIPS_history = sprintf('Z-score trials within subject column-wise\n%s\nRegressed out WM+CSF\nexplained %3.0f%% of variance on average', descrip, 100*mean(vexplained) );
ST_cleaned.SIIPS_history = textwrap({ST_cleaned.SIIPS_history}, 80);

figure; imagesc(mydata); colorbar
title('SIIPS')

% -----------------------------------------------------

%% PDM1
% -----------------------------------------------------
% Z-score within-subject

mydata = cellfun(@windsorize_matrix_columnwise, ST_data_PDM1_response, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cellfun(@zscore, mydata, 'UniformOutput', false); % Z-score within subject
mydata = cat(1, mydata{:});

% Regress out white-matter and CSF averages (doesn't make much difference)
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

ST_cleaned.pdm1 = mydata;
ST_cleaned.pdm1_history = sprintf('Z-score trials within subject column-wise\n%s\nRegressed out WM+CSF\nexplained %3.0f%% of variance on average', descrip, 100*mean(vexplained) );
ST_cleaned.pdm1_history = textwrap({ST_cleaned.pdm1_history}, 80);

figure; imagesc(mydata); colorbar
title('PDM1')

%% Region averages
% -----------------------------------------------------
% Z-score within-subject

mydata = cellfun(@windsorize_matrix_columnwise, ST_data_region_averages, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
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

ST_cleaned.big_regions = mydata;
ST_cleaned.big_regions_history = sprintf('Z-score trials within subject column-wise\n%s\nRegressed out WM+CSF\nexplained %3.0f%% of variance on average', descrip, 100*mean(vexplained) );
ST_cleaned.big_regions_history = textwrap({ST_cleaned.big_regions_history}, 80);

figure; imagesc(mydata); colorbar
title('Regions')

%% Region averages - fine-grained
% -----------------------------------------------------
% Z-score within-subject

mydata = cellfun(@windsorize_matrix_columnwise, ST_data_painpathways_finegrained, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
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

ST_cleaned.fine_regions = mydata;
ST_cleaned.fine_regions_history = sprintf('Z-score trials within subject column-wise\n%s\nRegressed out WM+CSF\nexplained %3.0f%% of variance on average', descrip, 100*mean(vexplained) );
ST_cleaned.fine_regions_history = textwrap({ST_cleaned.fine_regions_history}, 80);

figure; imagesc(mydata); colorbar
title('Fine-grained Regions')

%% Pain and temperature

% Temperature is not in C for all studies (sometimes coded by level)
% Pain ratings are not on the same scale - z-score within-person for our purposes here.

% Note: there may be sound data from BMRK5 data in here, and at least one
% study is missing temps, so this may need more work to fully analyze
% relationships with pain, etc. But it should suffice for preliminary
% identification of pathways.

% some missing values - impute mean subject- and column-wise
mydata = ST_data_pain_temp;
for i = 1:length(mydata)
    
    % work on each column
    wh_nan = isnan(mydata{i}(:, 1)); % pain is nan
    temps = mydata{i}(:, 2);                        % could customize for temperatures, but not done yet...
    
    mydata{i}(wh_nan, 1) = nanmean(mydata{i}(:, 1)); % replace with overall mean
end

% windsorize and Z-score - pain only
mydata = cellfun(@windsorize_matrix_columnwise, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cellfun(@zscore, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cat(1, mydata{:});

ST_cleaned.pain_rating = mydata(:, 1);

%% Figure out why pain and brain trials don't match and fix it

% Count pain trials per subjects
mydata = ST_data_pain_temp;
sz = cellfun(@size, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
sz = cat(1, sz{:});
n_pain_trials = sz(:, 1);

% Count pain trials per subjects
mydata = ST_data_painpathways_finegrained;
sz = cellfun(@size, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
sz = cat(1, sz{:});
n_brain_trials = sz(:, 1);

[n_pain_trials n_brain_trials]

% one subject, Stephan's placebo study.
wh = find(n_pain_trials ~= n_brain_trials)
Subject_Table.Studyname(wh), Subject_Table.Matfile{wh}
mymatfile = Subject_Table.Matfile{wh}
tmp = load(mymatfile); size(tmp.dat.dat) % 46 trials.
mypain = ST_data_pain_temp{wh}; mypain = mypain(:, 1);
sum(~isnan(mypain)) % 46 non-NaN trials.
whtrials = ~isnan(mypain);

ST_data_pain_temp{wh} = ST_data_pain_temp{wh}(whtrials, :);

% Re-do the pain processing

mydata = ST_data_pain_temp;
for i = 1:length(mydata)
    
    % work on each column
    wh_nan = isnan(mydata{i}(:, 1)); % pain is nan
    temps = mydata{i}(:, 2);                        % could customize for temperatures, but not done yet...
    
    mydata{i}(wh_nan, 1) = nanmean(mydata{i}(:, 1)); % replace with overall mean
end

% windsorize and Z-score - pain only
mydata = cellfun(@windsorize_matrix_columnwise, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cellfun(@zscore, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cat(1, mydata{:});

ST_cleaned.pain_rating = mydata(:, 1);

%% Get subject ID and study ID column
sz = cellfun(@size, ST_data_painpathways_finegrained, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
sz = cat(1, sz{:});
ntrials = sz(:, 1);  % number of trials per person
Subject_Table.ntrials = ntrials;

sindx = {};
studyindx = {};

for i = 1:length(ntrials)
    
    sindx{i} = i * ones(ntrials(i, 1), 1);
    
    studyindx{i} = Subject_Table.Study(i) * ones(ntrials(i, 1), 1);
    
end

subj_idx = cat(1, sindx{:});

studyindx = cat(1, studyindx{:});

ST_cleaned.subj_idx = subj_idx;

ST_cleaned.studyindx = studyindx;


%% NOW do temperature, and make sure we have uniform scale
% adjust study-by-study as needed
% examine coding of temperature
% discovered that we will need to get rid of pain/sound sound trials in bmrk5

create_figure('temp by pain');

for i = 1:max(Subject_Table.Study)
    
    allstudypaintemp = cat(1, ST_data_pain_temp{Subject_Table.Study==i});
    
    cla
    scatter(allstudypaintemp(:, 2) + .25 * (rand(length(allstudypaintemp), 1) - .5), allstudypaintemp(:, 1));
    refline; axis auto;
    nancorr(allstudypaintemp)
    drawnow
    
    %a = input(sprintf('Study %3.0f %s\n', i, studynames{i}));
end

% coded 1 2 3 : Study 3 BMRK4. was 46, 47, 48. replace.
% coded -1 0 1: Study 4. IE:  was 46, 47, 48
% coded 1 2 : Study 12. REMI. was individually calibrated. methods says 4 temps used...
% half the trials are NaN: study 7. bmrk5_painsound.  are these sound trials?? if so get rid of them.
% missing: study 11 placebo_value_stephan
% constant: Study 9: RCP.

%%
% wh_subj = Subject_Table.Subject(Subject_Table.Study == 3); % subject id number

temps = cat(1, ST_data_pain_temp{:});
temps = double(temps(:, 2));

temps(ST_cleaned.studyindx == 3 & temps == 1) = 46;
temps(ST_cleaned.studyindx == 3 & temps == 2) = 47;
temps(ST_cleaned.studyindx == 3 & temps == 3) = 48;

temps(ST_cleaned.studyindx == 4 & temps == -1) = 46;
temps(ST_cleaned.studyindx == 4 & temps == 0) = 47;
temps(ST_cleaned.studyindx == 4 & temps == 1) = 48;

% Relative temperature may be most helpful, as durations and thermode sizes also
% differ across studies.
% Get relative temperature:

temps_rel = NaN .* ones(size(temps));

for i = 1:length(Subject_Table.Subject)
    
    wh_trials = ST_cleaned.subj_idx == i;
    
    mytemps = temps(wh_trials);
    
    mytemps(~isnan(mytemps)) = zscore(mytemps(~isnan(mytemps)));
    
    %mytemps = rankdata(mytemps) ./ length(mytemps);
    
    temps_rel(wh_trials) = mytemps;
end

ST_cleaned.studynames = studynames;

ST_cleaned.rel_temp = temps_rel;

% study 12, remi: nan-out absolute temp
temps(ST_cleaned.studyindx == 12 & temps == 1) = NaN; % individually coded...
temps(ST_cleaned.studyindx == 12 & temps == 2) = NaN;

ST_cleaned.abs_temp = temps;

% study 11, placebo stephan: all temps are constant
wh = ST_cleaned.studyindx == 11;
ST_cleaned.abs_temp(wh) = 46.4;  
ST_cleaned.rel_temp(wh) = 0;
%% Get rid of sound trials in bmrk5_painsound

%%%wasremoved = studyindx == 7 & isnan(temps)

wh = ST_cleaned.studyindx == 7 & isnan(ST_cleaned.abs_temp);

N = fieldnames(ST_cleaned);
for i = 1:length(N)
    
    tmp = ST_cleaned.(N{i});
    
    if size(tmp, 1) == 27174 % original number of trials; we cut down to 25014 trials after removing sound trials
        
        tmp(wh, :) = [];
        
        ST_cleaned.(N{i}) = tmp;
    end
    
end

ST_cleaned.wasremoved = wh;
ST_cleaned.wasremoved_descrip = 'Sound trials removed from BMRK5; still in ST_* variables';

%% Plot pain x temp correlations/scatterplots


pain_temp_r = [];

create_figure('temp by pain', 4, 4);

for i = 1:max(Subject_Table.Study)
    
    wh = ST_cleaned.studyindx == i & ~isnan(ST_cleaned.rel_temp) & ~isnan(ST_cleaned.pain_rating);
    
    subplot(4, 4, i)
    scatter(ST_cleaned.rel_temp(wh) + .25 * (rand(length(ST_cleaned.rel_temp(wh)), 1) - .5), ST_cleaned.pain_rating(wh));
    refline; axis auto;
    drawnow
    
    if sum(wh)
        myr = corr([ST_cleaned.rel_temp(wh) ST_cleaned.pain_rating(wh)]);
        pain_temp_r(i, 1) = myr(1, 2);
    else
        pain_temp_r(i, 1) = NaN;
    end
    
    title(sprintf('Study %i %s', i, ST_cleaned.studynames{i}));
    
    %a = input(sprintf('Study %3.0f %s\n', i, studynames{i}));
end

Studynames = studynames';
pain_temp_r_table = table(Studynames, pain_temp_r);
disp(pain_temp_r_table);

saveas(gcf, 'figures/Pain_temperature_relationships_by_study.png');

%% Add gray/white/CSF

%ST_cleaned.gray_white_csf = cat(1, ST_data_gray_white_csf{:});
%ST_cleaned.gray_white_csf = ST_cleaned.gray_white_csf(~ST_cleaned.wasremoved, :);

mydata = ST_data_gray_white_csf;
mydata = cellfun(@windsorize_matrix_columnwise, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cellfun(@zscore, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cat(1, mydata{:});
mydata = mydata(~ST_cleaned.wasremoved, :);

ST_cleaned.gray_white_csf = mydata;

mydata = ST_data_gray_white_csf_components;
mydata = cellfun(@windsorize_matrix_columnwise, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cellfun(@zscore, mydata, 'UniformOutput', false); % Windsorize to 3 SD before Z-scoring
mydata = cat(1, mydata{:});
mydata = mydata(~ST_cleaned.wasremoved, :);

ST_cleaned.gray_white_csf_components = mydata;

%% Add labels

ST_cleaned.nps_labels = {pain_regions_nps.shorttitle};
ST_cleaned.siips_labels = {pain_regions_siips.shorttitle};

ST_cleaned.pdm1_labels = {pain_regions_pdm1.shorttitle};
ST_cleaned.fine_labels = pain_pathways_finegrained.labels;
ST_cleaned.big_labels = pain_pathways.labels;
ST_cleaned.graywhitecsf_labels = {'Gray_mean' 'White_mean' 'CSF_mean'};

ST_cleaned.graywhitecsfcomponents_labels = {'Gray1' 'Gray2' 'Gray3' 'Gray4' 'Gray5' 'White1' 'White2' 'White3' 'White4' 'White5' 'CSF1' 'CSF2' 'CSF3' 'CSF4' 'CSF5'};



%% Save

save data/pain_pathways_single_trial_data_2019 -append ST_cleaned Subject_Table





