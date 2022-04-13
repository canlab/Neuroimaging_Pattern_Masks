% Given an ROI defined by a mask, calculate seed co-activation map using the Neurosynth database
%
% :Usage:
% ::
%
%     [region_obj, stats] = neurosynth_interregion_coactivation(neurosynth_data, roi_mask);
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
%        (1) an atlas class object defining a series of regions
%
%        (2) A region-class object with one or more regions (elements)
%
%         - The atlas/mask object does not have to be in the same space/voxel size as the neurosynth_data object.
%         - It will be resampled to the space of neurosynth_data
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
%   **noplot:**
%        Not implemented yet.
%
% :Outputs:
% ::
%   **region_obj:**
%        region object with activation vectors for each Neurosynth contrast
%        in region_obj(i).dat
%
%        - Activation defined as any activation within ROI for each contrast
%
%   **stats:**
%        Structure with several fields:
%
%                    .r     Kendall's Tau b co-activation matrix                                                      
%                    .t     T-score matrix 
%                    .p,    P-value matrix 
%                    .fdrp  FDR-corrected q values
%                    .fdrthresh  P-value threshold for q < 0.05 FDR-corrected
%                    .sig   Significant positive co-activation with FDR correction
%            region_table:  Table object with region names and coordinates
%      distance_matrix_mm:  Distances between region centers in mm
%          distance_table:  Table object with pairwise correlation (phi,
%                           distance, and adjusted (residualized) distances
%     distance_fit_method:  Fit method for determining distance effect to residualize; default 'exp2'
%     r_distance_adjusted:  Adjusted (residualized) correlations among regions
%
% Notes on inferential stats:
% - Pearson's r gives the same point estimates as Kendall's tau b or Phi
%      coefficients, but inferential stats are different across these. 
% - Phi is dramatically faster (about 1000x) and point estimates are identical to Kendall's Tau b.
% - Phi inferential stats are more overconservative than Kendall's Tau b - in
%      one test, by an average of t(tau) - t(phi) = 7.48 (4/2022 Pain Pathways
%      intercorrelations), but t-values are 17.6 on average, so power is not an
%      issue. Phi seems appropriate and very fast to calculate.
%
% :Examples:
% ::
% ------------------------------------------------------------------------
% Get co-activation matrix for voxels in "pain pathways"
% ------------------------------------------------------------------------
%
% load neurosynth_data_obj.mat  % load the prepared data object
%
% atlas_obj = load_atlas('painpathways');
%
% [region_obj, stats] = neurosynth_interregion_coactivation(neurosynth_data, atlas_obj);
%
%
% ------------------------------------------------------------------------
% Get co-activation matrix for voxels in entire canlab 2018 atlas
% ------------------------------------------------------------------------
%
% load neurosynth_data_obj.mat  % load the prepared data object
%
% atlas_obj = load_atlas('canlab2018_2mm');
%
% [region_obj, stats] = neurosynth_interregion_coactivation(neurosynth_data, atlas_obj);
%
% This was used to create the file:
% save neurosynth_interregion_coactivation_canlab2018_2mm region_obj stats
%
% ------------------------------------------------------------------------
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

function [region_obj, stats] = neurosynth_interregion_coactivation(neurosynth_data, roi_mask, varargin)

doplot = true;

if any(strcmp(varargin, 'noplots')) || any(strcmp(varargin, 'noplot')), doplot = false; end


% ---------------------------------------------------------------------
% Conversion, error checking, and resampling
% ---------------------------------------------------------------------

% Convert to mask if region object is entered
if isa(roi_mask, 'region')

    roi_mask = region2atlas(roi_mask);

end

% Error checking

roi_mask = check_properties(roi_mask);

% Resampling
roi_mask = resample_space(roi_mask, neurosynth_data);

% atlas object should be int32 type, so no interpolation issues possible

% ---------------------------------------------------------------------
% Get roi_activation for each region
% ---------------------------------------------------------------------

disp('Getting activation vector across Neurosynth contrasts for each region')

u = unique(roi_mask.dat);
u(u == 0) = [];             % 0 is treated as empty/excluded here.

c = size(neurosynth_data.dat, 2);   % num contrasts
k = length(u);                      % num regions

roi_activation = false(c, k);

for i = 1:length(u)

    % find in-mask voxels
    % basically finding "true", but below is more robust if mask contains multiple integers or close-to-zero non-zeros
    wh = roi_mask.dat == u(i);

    % get data for each voxel in the ROI
    roi_data = neurosynth_data.dat(wh, :);

    % calculate whether each study activated within the ROI
    roi_activation(:, i) = any(roi_data, 1)';

end

% ---------------------------------------------------------------------
% calculate correlation matrix
% ---------------------------------------------------------------------

disp('Getting correlation matrix - Phi correlation')

% Note: Point estimates for Pearson's r, Kendall's Tau b, and Phi are
% numerically identical for binary data ({0, 1}).  The inferential
% statistics will differ among them.
%
% Phi is dramatically faster (about 1000x) and point estimates are identical to Kendall's Tau b.
% Phi inferential stats are more overconservative than Kendall's Tau b - in
% one test, by an average of t(tau) - t(phi) = 7.48 (4/2022 Pain Pathways
% intercorrelations), but t-values are 17.6 on average, so power is not an
% issue. Phi seems appropriate and very fast to calculate.

tic
[stats.r, stats.t, stats.p, stats.fdrp, stats.fdrthresh] = correlation('phi', roi_activation);
toc

% sig: positive co-activation only
stats.sig = stats.p < stats.fdrthresh & stats.t > 0;

% Tau calculates in 337 sec vs. 0.37 sec or so for Phi. 
% tic
% [stats.r, stats.t, stats.p, stats.fdrp, stats.fdrthresh] = correlation('taub', roi_activation);
% toc

% ---------------------------------------------------------------------
% Output
% ---------------------------------------------------------------------

% Save region object with data appended
disp('Appending ROI activation vectors to region object r(i).dat field')

region_obj = atlas2region(roi_mask);

for i = 1:k

    region_obj(i).dat = roi_activation(:, k);

end

% ---------------------------------------------------------------------
% Correction for distance
% ---------------------------------------------------------------------

region_mm_center = cat(1, region_obj.mm_center);

stats.region_table = table(roi_mask.labels', region_mm_center, 'VariableNames', {'Region' 'XYZmm_center'});

stats.distance_matrix_mm = squareform(pdist(region_mm_center));

stats.distance_table = table(squareform(stats.r - eye(size(stats.r)))', squareform(stats.distance_matrix_mm)', 'VariableNames', {'Phi_corr' 'Dist_in_mm'});

x = table2array(stats.distance_table);

stats.distance_fit_method = 'exp2';                             % bi-exponential
curve = fit(x(:, 2), x(:, 1), stats.distance_fit_method);       % x is distance, y is coactivation

[curve, goodness, output] = fit(x(:, 2), x(:, 1), 'exp2');

stats.distance_table.adjusted = output.residuals;
stats.r_distance_adjusted = squareform(stats.distance_table.adjusted);

xlabel('Distance in mm'); 
ylabel('Co-activation')

% ---------------------------------------------------------------------
% Plot
% ---------------------------------------------------------------------

if doplot

    create_figure('Neurosynth co-activation', 1, 3)
    plot_correlation_matrix(stats, 'var_names', roi_mask.labels, 'nofigure');
    title('Neurosynth co-activation')
    axis image; axis equal

    subplot(1, 3, 2)
    plot(curve, x(:, 2), x(:, 1));

    xlabel('Distance between regions (mm)');
    ylabel('Co-activation (phi)');
    title('Co-activation as f(distance)')

    subplot(1, 3, 3)
    % create temporary copy to plot adjusted
    stats_tmp = stats;
    stats_tmp.r = stats.r_distance_adjusted;

    
%     clim = max(abs(stats.r_distance_adjusted(:)));

    plot_correlation_matrix(stats_tmp, 'var_names', roi_mask.labels, 'nofigure'); % , 'colorlimit', [-clim clim]);
    title('Co-act, Distance-corrected')
    axis image; axis equal
end



end % main function


% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% SUBFUNCTIONS
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
