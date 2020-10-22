function yhat = classify_somatovisceral_pain(dat)
%classify_somatovisceral_pain: Apply network-based classification model
%that discriminates between somatic and visceral pain from Van Oudenhove,
%Kragel et al. 2020 Nature Communications.
%
% This function makes a specific kind of ROC curve plot, based on input
% values along a continuous distribution and a binary outcome variable
% (logical)
%
% :Usage:
% ::
%
%   yhat = classify_somatovisceral_pain(dat)
%
% :Inputs:
%   dat - fmri data object to be classified
%
% :Optional Inputs:
%   -none for now-
%
% :Outputs:
%
%   Yhat - a vector with one continous estimate of the probability an image
%   is somatic (positive values) or visceral pain (values near zero)
%
%   Uses the model parameters saved in Visceral_vs_Somatic_betas_Yeo_Networks.mat
%
% :Examples:
% ::
%
%    dat = load_image_set('bmrk3');
%    yhat = classify_somatovisceral_pain(dat);
% ..
%    Philip Kragel Oct 2020
%    See Notes in text for more programming details.
%

%    Notes:
%
%    10/22/2020: Phil Kragel - initial version from repo code
%
% ..
load(which('Visceral_vs_Somatic_betas_Yeo_Networks.mat'),'b_avg')
rsfmri_stats=image_similarity_plot(dat,'average','mapset','bucknerlab','noplot');
r=[rsfmri_stats.r];
Z=fisherz(r)';
yhat = zeros(size(dat.dat,2),7);
for j=1:7
    yhat(:,j)=glmval(b_avg(j,:)',Z(:,j),'logit');
    
end
yhat=nanmean(yhat,2);