function mean_data = mean_data_cross_trials(meta, regional_data, study_canlab_dataset)
% meta and regional_data are the outputs of function preparation_extract_data
% remove high vif trials
% regress global mean, WM and CSF effects from regional average activity
% categorize by temp,calculate mean activity

% ----------------Rating rescale-----------------------------
Pain = single_trial_retrieve_data_all_studies(study_canlab_dataset, 'ratings');
Pain = single_trial_rescale_multistudy_data_by_mad(Pain);
idx = [ 9, 11, 13]; % index of removed studies
Pain.event_by_study_madscaled(idx) = [];
Pain.event_by_study(idx) = [];

Nstudies = length(regional_data);
study_idx = 1:Nstudies;

% remove placebo/drug/high-vif trials 
%--------------------------------------------------------------------------
pd_idx = cellfun(@isempty, meta.placebo_drug_trial);
pd_idx = find(pd_idx == 0);

for i = 1:Nstudies % remove high-vif trials
    
    Nsubj = length(meta.meta_data{i});
           
    for j = 1:Nsubj
        
        Hvif_idx = meta.meta_data{i}{j}(:,6) == 0; % find ok trials = 0
        regional_data{i}{j}(Hvif_idx, :) = [];
        meta.meta_data{i}{j}(Hvif_idx,:) = []; 
%         meta.meta_data{i}{j}(:, 1) = Pain.event_by_study_madscaled{i}{j};
        meta.meta_data{i}{j}(:, 1) = Pain.event_by_study{i}{j};
        
        if ismember(study_idx(i), pd_idx) % if drug/placebo exist, remove these trials  
        
           meta.placebo_drug_trial{i}{j}(Hvif_idx,:) = []; 
           Drug_idx = meta.placebo_drug_trial{i}{j} == 1; 
           regional_data{i}{j}(Drug_idx, :) = [];
%            meta.meta_data{i}{j}(Drug_idx,:) = []; 
           meta.meta_data{i}{j}(Drug_idx,:) = [];
        
        end 
               
    end

end
   
% calculate the mean activity cross trials
%--------------------------------------------------------------------------

for i = 1:Nstudies
    
    Nsubj = length(meta.meta_data{i});
    
    for j = 1:Nsubj
        
        temp_tmp = meta.meta_data{i}{j}(:,2); % temperature
        temp_deg = unique(temp_tmp); % degree of temperature level
        temp_deg(isnan(temp_deg)) = [];
        temp_level{i}{j} = temp_deg;
        
        for t = 1:length(temp_deg)
            
            idx = find(temp_tmp == temp_deg(t)); % index of trials under temperature level t
            A = meta.meta_data{i}{j}(idx, 1); 
            A(isnan(A)) = []; % remove nan trials
            ratings_mean{i}{j}(t,:) = mean(A); % average across trials
                     
%             activity_mean{i}{j}(t,:) = mean(data_here{i}{j}(idx, :));  
            activity_ori_mean{i}{j}(t,:) = mean(regional_data{i}{j}(idx, :)); % original averaged activity
            RegressData{i}{j}(t,:) = mean(meta.meta_data{i}{j}(idx, [3 4 5])); % signials needs to regress (gloable GM, WM and CSF signal)
                        
        end   
                
    end
    
end
mean_data.rating = ratings_mean;
mean_data.temp_level = temp_level;
mean_data.activity_ori_mean = activity_ori_mean;
mean_data.regress_var = RegressData;

% linear regression
% -------------------------------------------------------------------------
for i = 1:Nstudies
    
    EventDataSeq = cat(1,mean_data.activity_ori_mean{i}{:}); % 
    Regressor = cat(1,mean_data.regress_var{i}{:});  % regress global mean, WM, CSF signals
    
        for k = 1:size(EventDataSeq, 2) % number of regions
            
            X = [ones(size(Regressor,1),1) Regressor];
            X1 = (X'*X)\X'; X1 = X1 * EventDataSeq(:, k);
            data_r(:, k) = EventDataSeq(:, k) - X*X1;
%            [b, bint, data_r(:,k)] = regress(EventDataSeq(:, k),X);           
                        
        end  
        
        for j = 1:length(mean_data.activity_ori_mean{i}) % Nsubj
            
            Nmean = size(mean_data.activity_ori_mean{i}{j}, 1); 
            activity_mean_r{i}{j} = data_r(1:Nmean, :); % averaged activity with global mean, WM, CSF signals are regressed
            data_r(1:Nmean, :) = [];
            
        end
        
        clear data_r;
    
end

mean_data.activity_mean_r = activity_mean_r;

% mad rescale
% -------------------------------------------------------------------------
for i = 1:Nstudies
    
    % activity    
    A = cat(1,mean_data.activity_mean_r{i}{:}); 
    B = cat(1,mean_data.activity_ori_mean{i}{:});
  
    for j = 1:size(A, 2) % Nregions
        
        Mad_r = mad(A(:, j));
        Mad_ori = mad(B(:, j));
        
        for k = 1:length(mean_data.activity_mean_r{i}) %Nsubj
    
            activity_mean_mad{i}{k}(:,j) = mean_data.activity_mean_r{i}{k}(:,j) ./ Mad_r;     
            activity_ori_mean_mad{i}{k}(:,j) = mean_data.activity_ori_mean{i}{k}(:,j) ./ Mad_ori;
            
        end
        
    end 
    
end
mean_data.activity_mean_r_mad = activity_mean_mad;
mean_data.activity_ori_mean_mad = activity_ori_mean_mad;
mean_data.study_name = meta.study_name;
