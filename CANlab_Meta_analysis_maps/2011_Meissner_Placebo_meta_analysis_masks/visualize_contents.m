function visualize_contents()
% visualize_contents  Render Meissner 2011 placebo convergence masks.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = {
    'Meissner2011_PlaceboInc_10mm', 'placebo_increases_at_least_3_studies_in_10_mm.hdr'
    'Meissner2011_PlaceboInc_15mm', 'placebo_increases_at_least_3_studies_in_15_mm.hdr'
    'Meissner2011_PlaceboDec_10mm', 'placebo_decreases_at_least_3_studies_in_10_mm.hdr'
};

canlab_render_patterns(imgs);

end
