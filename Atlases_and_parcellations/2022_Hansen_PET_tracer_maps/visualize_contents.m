function visualize_contents()
% visualize_contents  Render Hansen 2022 PET-tracer neurotransmitter
% maps into png_images/. Tracers are continuous fmri_data patterns
% (not a labelled atlas), so we use canlab_render_patterns. See
% contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Pull the combined / z-scored / averaged set used by
% load_image_set('hansen22') and write one PNG per target.
[obj, target_names] = load_image_set('hansen22');

savedir = fullfile(this_dir, 'png_images');
if ~exist(savedir, 'dir'), mkdir(savedir); end

for k = 1:size(obj.dat, 2)
    label = matlab.lang.makeValidName(['hansen22_' target_names{k}]);
    sub   = get_wh_image(obj, k);

    try
        fig = create_figure([label '_surface']);
        axis off;
        try
            surface(sub, 'foursurfaces_hcp', 'nolegend', 'nofigure');
        catch
            clf(fig); surface(sub, 'nolegend', 'nofigure');
        end
        drawnow;
        saveas(fig, fullfile(savedir, [label '_surface.png']));
    catch ME
        warning('hansen_pet:surface', '%s: %s', label, ME.message);
    end

    try
        fig = create_figure([label '_montage']);
        axis off;
        montage(sub, 'nolegend');
        drawnow;
        saveas(fig, fullfile(savedir, [label '_montage.png']));
    catch ME
        warning('hansen_pet:montage', '%s: %s', label, ME.message);
    end
end

fprintf('\nDone. PNGs written to: %s\n', savedir);

end
