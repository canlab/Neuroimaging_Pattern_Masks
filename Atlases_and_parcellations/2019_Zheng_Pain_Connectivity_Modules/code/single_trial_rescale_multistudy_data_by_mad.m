function VAR = single_trial_rescale_multistudy_data_by_mad(VAR)
%
% Rescale data for display/aggregation/analysis
% Both pain reports and brain data are on different scales across studies
% Rescale all data by single, study-level scaling factor based on median
% abs. deviation across observations.
% 
% Example:
% PAIN = single_trial_rescale_multistudy_data_by_mad(PAIN)
%
% dACC NPS subregion
% dACC = single_trial_retrieve_data_all_studies(study_canlab_dataset, 'dACC');
% dACC = single_trial_rescale_multistudy_data_by_mad(dACC);

% Nstudies = length(VAR.event_by_study);
Nstudies = length(VAR.event_by_study);

for i = 1:Nstudies
    
    myevents = VAR.event_by_study{i};
    
    eventscat = myevents;
    
    study_mad = mad(double(cat(1, eventscat{:})));
    
    % for each subject
    
    for j = 1:length(myevents)
        
        myevents{j} = myevents{j} ./ study_mad;
        
    end % subject loop
    
    VAR.event_by_study_madscaled{i} = myevents;
    
end % study loop

tmp = cat(2, VAR.event_by_study_madscaled{:});
VAR.events_cat_madscaled = cat(1, tmp{:});

end % function

