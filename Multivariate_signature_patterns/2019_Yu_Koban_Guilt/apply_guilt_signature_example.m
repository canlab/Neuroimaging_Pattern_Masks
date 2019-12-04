% Applying a signature to new test data
% The guilt weight map is Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii.gz in the repository.

% If you have our "CanlabCore" repository from Github on your Matlab path, 
% an easy way to apply the signature pattern to a new dataset is the "apply_mask" 
% method for fmri_data objects.  For example, this code will load a sample dataset 
% and apply the signature.
% It returns a vector of one pattern response per input image.

test_images = load_image_set('emotionreg');                 % Load 30 emotion regulation-related contrast maps from Wager et al. 2008, Neuron

guilt = load_image_set('guilt');  % Load the Guilt Behavior map as an fmri_data object

guilt_response = apply_mask(test_images, guilt, 'pattern_expression'); % Apply the weight map

figure; plot(guilt_response, 'o'); title('Guilt response scores for each participant');
