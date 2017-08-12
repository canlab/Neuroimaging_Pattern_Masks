function [siips_pos_obj, siips_neg_obj, siips_pos_regions, siips_neg_regions, names_pos, names_neg, r] = load_siips_subregions()
%
% [siips_pos_obj, siips_neg_obj, siips_pos, siips_neg, r] = load_siips_subregions()
%
% Load the subregions from the SIIPS1 pattern image.
% SIIPS1 is described in Woo et al. 2017, Nature Communications.
%
% siips_pos_obj, siips_neg_obj : fmri_data objects containing local pattern weights for positive
% and negative subregions of the SIIPS1 pattern.  region(siips_pos_obj) returns a region object with a vector
% of local patterns.
%
% siips_pos_regions, siips_neg_regions : region objects containing local
% pattern weights (in siips_pos_regions(x).all_data, etc.) and region names
% (siips_pos_regions.shorttitle)
%
% names_pos, names_neg : names for the subregions in the structures above.
% 
% r : a region object with all names for positive and negative subregions,
% in r(x).shorttitle


% get regions
% -------------------------------------------------------------------------

%siips_parcel_object = fmri_data(which('nonnoc_v11_4_subcluster_maps_fdr05_unique_wttest.nii'), 'noverbose');

% Using the patterns works because we have defined regions with contiguous
% values in this case. Then we can return local patterns in region object.
siips_pattern_object = fmri_data(which('nonnoc_v11_4_subcluster_maps_fdr05_pattern_wttest.nii'), 'noverbose');

% r(x).all_data contains local weight patterns.
r = region(siips_pattern_object, siips_pattern_object, 'noverbose');

% get names
% -------------------------------------------------------------------------

names = load(which('nonnoc_v11_4_subcluster_maps_fdr05_44_cluster_names.mat'));

for i = 1:length(r)
    
    r(i).shorttitle = deblank(names.cl44names(i, :)); 
    
end


% separate positive and negative regions
% -------------------------------------------------------------------------

whpos = names.clsign > 0;
whneg = names.clsign < 0;

siips_pos_regions = r(whpos);
siips_neg_regions = r(whneg);

names_pos = names.cl44names(whpos, :);
names_neg = names.cl44names(whneg, :);

for i = 1:length(siips_pos_regions)

    siips_pos_regions(i).shorttitle = deblank(names_pos(i, :)); 
    
    siips_pos_regions(i).descrip1 = sprintf('SIIPS1, index %d among pain-positive subregions', i);
end

for i = 1:length(siips_neg_regions)

    siips_neg_regions(i).shorttitle = deblank(names_neg(i, :)); 
    
    siips_neg_regions(i).descrip1 = sprintf('SIIPS1, index %d among pain-negative subregions', i);
end

% Get image objects for pos and neg with local patterns
% -------------------------------------------------------------------------
tmpmask = region2fmri_data(siips_pos_regions, siips_pattern_object); 

% this still returns parcel indices, until we add functionality to
% region2fmri_data. so we need to mask the original pattern.

siips_pos_obj = apply_mask(siips_pattern_object, tmpmask);

tmpmask = region2fmri_data(siips_neg_regions, siips_pattern_object); 

% this still returns parcel indices, until we add functionality to
% region2fmri_data. so we need to mask the original pattern.

siips_neg_obj = apply_mask(siips_pattern_object, tmpmask);

% reparse contiguous voxels because extract_roi_averages needs this
siips_pos_obj = reparse_contiguous(siips_pos_obj, 'nonempty');
siips_neg_obj = reparse_contiguous(siips_neg_obj, 'nonempty');

end % main function

