function visualize_contents()
% visualize_contents  Render placebo-prediction patterns into png_images/.
% Wager et al. 2011 J Neurosci. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Make subfolders' .img files findable.
addpath(genpath(this_dir));

imgs = {
    'PlaceboPredict_Anticipation',     'PlaceboPredict_Anticipation.img'
    'PlaceboPredict_PainPeriod',       'PlaceboPredict_PainPeriod.img'
    'PlaceboBrainPredict_Anticipation','PlaceboBrainPredict_Anticipation.img'
    };

canlab_render_patterns(imgs);

end
