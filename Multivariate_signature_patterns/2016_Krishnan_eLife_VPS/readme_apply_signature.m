% Applying a signature to new test data
% The VPS weight map is Rating_Weights_LOSO_2.nii in the repository.

% If you have our ?CanlabCore? repository from Github on your Matlab path, 
% an easy way to apply the signature pattern to a new dataset is the ?apply_mask? 
% method for fmri_data objects.  For example, this code will load a sample dataset 
% and apply the 'Vicarious Pain Signature' (VPS). 
% It returns a vector of one VPS response per input image.

test_images = load_image_set('emotionreg');                 % Load 30 emotion regulation-related contrast maps from Wager et al. 2008, Neuron

pines_sig = fmri_data(which('bmrk4_VPS_unthresholded.nii'));  % Load the VPS map as an fmri_data object

pines_response = apply_mask(test_images, pines_sig, 'pattern_expression'); % Apply the weight map

figure; plot(pines_response, 'o'); title('VPS response scores for each participant');
