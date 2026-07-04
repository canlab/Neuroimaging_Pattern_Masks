function canlab_render_atlas(atlases, varargin)
% canlab_render_atlas  Render montage + isosurface PNGs for atlas objects.
%
% Sibling of canlab_render_patterns for atlas-type objects. Each entry of
% the input describes one atlas to render and saves
%   <label>_montage.png      axial+sagittal montage of atlas regions
%   <label>_isosurface.png   3-D isosurface coloured by parcel index
% into the chosen save directory.
%
% Usage:
%   canlab_render_atlas(atlases)
%   canlab_render_atlas(atlases, 'savedir', some_dir)
%
% Inputs:
%   atlases : N-by-2 cell array. Column 1 is the short label used for the
%             PNG basename. Column 2 is *either*:
%               - the keyword passed to load_atlas(), e.g. 'glasser_fsl6'
%               - a .mat filename on the path containing an `atlas_obj`
%                 variable (or `atlas` object directly)
%               - an absolute path to a .mat file
%
% Optional name/value:
%   'savedir' : output directory for PNGs (default: <pwd>/png_images).
%
% Requirements: CanlabCore on path, SPM12 on path.

p = inputParser;
p.addParameter('savedir', fullfile(pwd, 'png_images'), @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;

if ~exist(opts.savedir, 'dir'), mkdir(opts.savedir); end

for k = 1:size(atlases, 1)

    label = atlases{k, 1};
    spec  = atlases{k, 2};

    fprintf('\n=== %s ===\n  %s\n', label, spec);

    atlas_obj = [];
    try
        % First try as a load_atlas keyword.
        atlas_obj = load_atlas(spec);
    catch
        % Then try to resolve as a .mat file containing an atlas.
        fpath = which(spec);
        if isempty(fpath) && exist(spec, 'file')
            fpath = spec;
        end
        if isempty(fpath)
            warning('canlab_render_atlas:missing', ...
                'Could not resolve atlas spec %s (not a load_atlas keyword and not on path).', spec);
            continue
        end
        S = load(fpath);
        % Pick the first variable that looks like an atlas.
        fn = fieldnames(S);
        for j = 1:numel(fn)
            v = S.(fn{j});
            if isa(v, 'atlas')
                atlas_obj = v; break
            end
        end
        if isempty(atlas_obj)
            warning('canlab_render_atlas:notatlas', ...
                'No atlas-typed variable found in %s.', fpath);
            continue
        end
    end

    % --- Axial+sagittal montage of atlas regions ----------------------
    try
        r = atlas2region(atlas_obj);
        fig = figure('Name', [label '_montage']);
        o2 = canlab_results_fmridisplay([], 'multirow', 1);
        try, brighten(0.6); end %#ok<TRYNC>
        o2 = montage(r, o2, 'wh_montages', 1:2);
        try, scn_export_papersetup(600); end %#ok<TRYNC>
        drawnow;
        saveas(fig, fullfile(opts.savedir, [label '_montage.png']));
    catch ME
        warning('canlab_render_atlas:montage', ...
            '%s: montage failed (%s)', label, ME.message);
    end

    % --- 3-D isosurface ----------------------------------------------
    try
        fig = figure('Name', [label '_isosurface']);
        han = isosurface(atlas_obj);
        if ~isempty(han), set(han, 'FaceAlpha', 0.5); end
        view(135, 20);
        try, lightFollowView; end %#ok<TRYNC>
        try, lightRestoreSingle; end %#ok<TRYNC>
        axis off; drawnow;
        saveas(fig, fullfile(opts.savedir, [label '_isosurface.png']));
    catch ME
        warning('canlab_render_atlas:isosurface', ...
            '%s: isosurface failed (%s)', label, ME.message);
    end

end

fprintf('\nDone. PNGs written to: %s\n', opts.savedir);

end % function
