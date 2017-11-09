% Copyright 2016 Monica Rosenberg and Emily Finn

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% Rosenberg MD, Finn ES, Scheinost D, Papademetris X, Shen X, 
% Constable RT & Chun MM. (2016). A neuromarker of sustained
% attention from whole-brain functional connectivity. Nature Neuroscience 
% 19(1), 165-171.

% This code provides a framework for testing the sustained attention 
% connectome-based predictive model (sustained attention CPM) in a novel 
% dataset, as described in Rosenberg, Finn et al 2016 (see above for full 
% reference). The input 'validation_mats' is a pre-calculated 268x268xN matrix 
% containing all individual-subject connectivity matrices in the test set,
% where 268 is number of nodes in the Shen 2013 brain atlas (available 
% https://www.nitrc.org/frs/?group_id=51) and N = number of subjects. Each 
% element (i,j,k) in these matrices represents the correlation between the 
% BOLD timecourses of nodes i and j in subject k during a single fMRI 
% session. The second input ('validation_behav') is the Nx1 vector of 
% scores for the behavior of interest for all subjects.

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.

%% Define variables
clear; clc;

% Load network masks
load('saCPM.mat')

% Validation data
validation_mats  = ;                        % validation data (n_node x n_node x n_validation_sub symmetrical connectivity matrices)
validation_behav = ;                        % n_validation_sub x 1 vector of behavior
n_validation_sub = size(validation_mats,3); % total number of validation subjects
   
%% Test model

for k = 1:n_validation_sub
    % Get network strength for every subject
    high_attn_strength(k,1) = sum(sum(high_attention_mask.*validation_mats(:,:,k)));
    low_attn_strength(k,1)  = sum(sum(low_attention_mask.*validation_mats(:,:,k)));
    
    % Generate predicted d' (a measure of how subjects would hypothetically 
    % perform on the gradCPT, a test of sustained attention) using the
    % original model from Rosenberg, Finn et al. 2016
    glm_pred_orig(k,1) = robGLM_fit_orig(1) + robGLM_fit_orig(2)*high_attn_strength(k) + robGLM_fit_orig(3)*low_attn_strength(k);
 
    % Generate predicted d' using the updated model, which includes one
    % predictor (low-attention network strength subtracted from 
    % high-attention network strength)
    glm_pred_update(k,1) = robGLM_fit_update(1) + robGLM_fit_update(2)*sum(sum(saCPM_mask.*validation_mats(:,:,k)));

end

%% Evaluate performance

% Correlate predicted and observed behavior
[r_pos,        p_pos]        = corr(validation_behav, high_attn_strength);
[r_neg,        p_neg]        = corr(validation_behav, low_attn_strength);
[r_glm_orig,   p_glm_orig]   = corr(validation_behav, glm_pred_orig);
[r_glm_update, p_glm_update] = corr(validation_behav, glm_pred_update);

% Plot results
ax1 = subplot(2,2,1);
scatter(validation_behav, high_attn_strength); lsline
title(['High-attention r = ' num2str(round(r_pos*100)/100) ', p = ' num2str(p_pos)])
xlabel('Observed'); ylabel('Predicted')

ax2 = subplot(2,2,2);
scatter(validation_behav, low_attn_strength); lsline
title(['Low-attention r = ' num2str(round(r_neg*100)/100) ', p = ' num2str(p_neg)])
xlabel('Observed'); ylabel('Predicted')

ax3 = subplot(2,2,3);
scatter(validation_behav, glm_pred_orig); lsline
title(['Original glm r = ' num2str(round(r_glm_orig*100)/100) ', p = ' num2str(p_glm_orig)])
xlabel('Observed'); ylabel('Predicted')

ax4 = subplot(2,2,4);
scatter(validation_behav, glm_pred_update); lsline
title(['Updated glm r = ' num2str(round(r_glm_update*100)/100) ', p = ' num2str(p_glm_update)])
xlabel('Observed'); ylabel('Predicted')
