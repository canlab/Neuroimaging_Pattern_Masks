% Given an ROI defined by a mask, calculate seed co-activation map using the Neurosynth database
%
% :Usage:
% ::
%
%     [r_map, t_map] = neurosynth_seed_coactivation_map(neurosynth_data, region_mask);
%
% :Inputs:
% ::
%   **neurosynth_data:**
%        The Neurosynth dataset formatted as an fmri_data object
%        Included in this distribution in the file neurosynth_data.mat 
%
%         - Neurosynth dataset is a [1, 0] matrix stored in an fmri_data object
%           called neurosynth_data
%           This was created from the MKDA analysis, using "unweighted_study_data" 
%         - MC_Info.mat has variable MC_Setup, which contains the dataset
%         - MC_Setup.unweighted_study_data is a voxels x study contrasts matrix of 1 and 0 values
%           1 means activation, 0 means no activation for each contrast.
%
%   **region mask:**
%        Two possible input formats:
%        (1) Any image_vector, fmri_data, or atlas class object with a data
%        field (.dat) consisting of values in the set {1,0}
%
%        (2) A region-class object with one or more regions (elements)
%
%         - The mask object does not have to be in the same space/voxel size as the neurosynth_data object.
%         - It will be resampled to the space of neurosynth_data
%         - If you enter a region object with multiple regions, they will
%         be "flattened", so that the "seed" will include all voxels in any
%         region. To get separate maps for each region, pass in regions one
%         at a time, e.g., region_obj(1), region_obj(2), etc.
%
%         Example:
%         rois = load_atlas('canlab2018_2mm');
%         vpl = select_atlas_subset(rois, {'VPL'});
%
%        This funtion calculates "activation" of the seed region this way:
%         - Define study contrasts-level activation as any activation within the
%           region of interest.
%         - The MKDA setup has already smoothed the peak reported coordinates with
%           a kernel (4 mm in the latest as of 4/2022)
%         - ACTIVATION thus means that a study produced a coordinate within 4 mm of
%           any voxel in the region.
%
% :Optional Inputs:
% ::
%   **param1:**
%        none yet.
%
% :Outputs:
% ::
%   **r_map:**
%        fmri_data object with correlation map (all voxels) in space of neurosynth data
%
%   **t_map:**
%        statistic_image object with t-score map (seed region removed)
% 
%        - resampled to space of neurosynth data.  
%        - Thresholded at 1e-8 (standard GWAS correction) for convenience, by default.
%        - Can re-threshold like so: t_map = threshold(t_map, 1e-6, 'unc');
%
%
% :Examples:
% ::
% ------------------------------------------------------------------------
% Get co-activation map for voxels in bilateral VPL thalamus
% ------------------------------------------------------------------------
%
% load neurosynth_data_obj.mat  % load the prepared data object
% 
% rois = load_atlas('canlab2018_2mm');
% 
% vpl = select_atlas_subset(rois, {'VPL'});
% 
% [r_map, t_map] = neurosynth_seed_coactivation_map(neurosynth_data, vpl);
% 
% orthviews(t_map);
%
% ------------------------------------------------------------------------
% Get co-activation map for voxels in bilateral VPL or VPM thalamus
% ------------------------------------------------------------------------
%
% load neurosynth_data_obj.mat  % load the prepared data object
% 
% rois = load_atlas('canlab2018_2mm');
% 
% vplm = select_atlas_subset(rois, {'VPL' 'VPM'});
% 
% % Convert to region object to demonstrate this functionality
% r = atlas2region(vplm);
%
% [r_map, t_map] = neurosynth_seed_coactivation_map(neurosynth_data, r);
% 
% orthviews(t_map);
%
% :References:
%   Yarkoni, Poldrack, Van Essen, & Wager 2011, Nature Methods
%
% :See also:
%   correlation.m, Meta_cluster_tools
%

% ..
%    Programmers' notes:
%     by Tor Wager, 4/2022
%     see also:
%     correlation.m, Meta_cluster_tools
%     [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,dat,[volInfo])
%     [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,R.dat,R.volInfo);
%     Meta_cluster_tools('activation_table', DB, MC_Setup, cl(1), 'Modality2', 'PosNeg');
% ..
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2022 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%

function [r_map, t_map] = neurosynth_seed_coactivation_map(neurosynth_data, roi_mask)

% ---------------------------------------------------------------------
% Conversion, error checking, and resampling
% ---------------------------------------------------------------------

% Convert to mask if region object is entered
if isa(roi_mask, 'region')

    roi_mask = region2atlas(roi_mask);

    % "flatten" to make into a mask
    roi_mask.dat(roi_mask.dat > 0) = 1;

end

% Error checking

if ~all(roi_mask.dat == 0 | roi_mask.dat == 1, 1) %any(double(unique(roi_mask.dat)) - [0 1]')
    error('roi_mask does not seem to contain only [1 0] values');
end

roi_mask = resample_space(roi_mask, neurosynth_data);

% just in case of interpolation issues, round
% intention is for mask to be [0, 1] only
roi_mask.dat = round(roi_mask.dat);

% ---------------------------------------------------------------------
% Correlation map and inferential statistics
% ---------------------------------------------------------------------

% find in-mask voxels
% basically finding "true", but below is more robust if mask contains multiple integers or close-to-zero non-zeros
wh = roi_mask.dat > 100*eps; 

% get data for each voxel in the ROI
roi_data = neurosynth_data.dat(wh, :);

% calculate whether each study activated within the ROI
roi_activation = any(roi_data, 1)';

% calculate correlation
r = get_large_data_correlation(neurosynth_data, roi_activation);

r_map = neurosynth_data;
r_map.dat = r;
r_map.fullpath = [];
r_map.image_names = [];

% Get T and P-values

r2t = @(r, n) r .* sqrt((n - 2) ./ (1 - r.^2));

t2p = @(t, n) 2 .* (1 - tcdf(abs(t), n - 2));

n = length(roi_activation);  % number of contrasts in neurosynth (could reduce df to be conservative)

t = r2t(r, n);
pp = t2p(t, n);

% Create t statistic image object and threshold
% remove seed voxels from map; they are invalid (circular)

t_map = statistic_image('dat', t, 'type', 't', 'dfe', n - 2, 'p', pp, 'volInfo', neurosynth_data.volInfo);
t_map.dat(wh, :) = 0;
t_map.p(wh, :) = 1;
t_map.sig(wh, :) = 0;

t_map = threshold(t_map, 1e-8, 'unc');  % Correct with standard threshold for GWAS analysis

% t_map = threshold(t_map, 0.05, 'fdr');

 
end % main function


% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% SUBFUNCTIONS
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

function r = get_large_data_correlation(neurosynth_data, roi_activation)


meanfunc = @(x) sum(x) ./ size(x, 1);  % fast mean for binary data

meancenter = @(x) x - meanfunc(x); 

% This version takes an N x p matrix a and an N x v matrix b and returns
% a p x v matrix of correlations across the pairs. 

corr_matrix = @(a, b) (meancenter(a)' * meancenter(b) ./ (size(a, 1) - 1)) ./ (std(b)' * std(a))'; % Correlation of a with each column of b

% Too large to fit in memory. 
% Calculate piecewise
% A chunk size of about 500 is optimal based on numerical tests, 4/2022

v = size(neurosynth_data.dat, 1);
st = 1:500:v;                       % starting indices
en = [st(2:end) - 1 v];             % ending indices
k = length(st);                     % chunks

r = NaN .* zeros(v, 1);

tic;
fprintf('Seed correlation in %3.0f chunks: %03d', k, 1)

for i = 1:k
    
    rr = corr_matrix(roi_activation, neurosynth_data.dat(st(i) : en(i), :)');
    
    r(st(i):en(i)) = rr;
    
    fprintf('\b\b\b%03d', i)
    
end
fprintf('\n')
toc

% t = [];
% xx = [10:100:2000];
% offset = 50000;
% for k = xx
%     t1 = clock;
%     r = corr_matrix(roi_activation, neurosynth_data.dat(offset + 1:offset + k, :)');
%     t(end+1) = etime(clock, t1);
%     
% end
% 
% % figure; 
% hold on; 
% plot(xx, t ./ xx);
% xlabel('Chunk size')
% ylabel('time per var');
% 
% % hold on; plot(xx, t(1) .* [1:10])

end
