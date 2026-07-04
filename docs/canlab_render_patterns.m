function canlab_render_patterns(imgs, varargin)
% canlab_render_patterns  Render surface / montage / isosurface PNGs for a list of NIfTIs.
%
% Used by the per-folder visualize_contents.m scripts in
% Neuroimaging_Pattern_Masks. Standardises output so that every study
% folder has the same triplet of figures in png_images/.
%
% Usage:
%   canlab_render_patterns(imgs)
%   canlab_render_patterns(imgs, 'savedir', some_dir)
%   canlab_render_patterns(imgs, 'kinds', {'surface', 'montage'})  % skip isosurface
%   canlab_render_patterns(imgs, 'thresh', 'fdr05')                % atlas-style isosurface
%
% Inputs:
%   imgs  : N-by-2 cell array. Column 1 is the short label used for the
%           PNG basename; column 2 is the NIfTI filename (.nii or .nii.gz)
%           resolved either on the MATLAB path or relative to the savedir.
%
% Optional name/value:
%   'savedir' : output directory for PNGs (default: <pwd>/png_images).
%   'kinds'   : cellstr of which renderings to produce; subset of
%               {'surface','montage','isosurface'}. Default: all three.
%   'thresh'  : optional threshold name passed to fmri_data threshold()
%               (default: '' = no threshold).
%
% Requirements: CanlabCore on path, SPM12 on path.

p = inputParser;
p.addParameter('savedir', fullfile(pwd, 'png_images'), @(x) ischar(x) || isstring(x));
p.addParameter('kinds', {'surface', 'montage', 'isosurface'}, @iscellstr);
p.addParameter('thresh', '', @ischar);
p.parse(varargin{:});
opts = p.Results;

if ~exist(opts.savedir, 'dir'), mkdir(opts.savedir); end

for k = 1:size(imgs, 1)

    label   = imgs{k, 1};
    relname = imgs{k, 2};

    fpath = which(relname);
    if isempty(fpath) || ~exist(fpath, 'file')
        % Look next to the caller for either the exact name or a .gz fallback.
        candidates = { fullfile(pwd, relname), fullfile(pwd, [relname '.gz']) };
        fpath = '';
        for c = 1:numel(candidates)
            if exist(candidates{c}, 'file')
                fpath = candidates{c};
                break
            end
        end
    end
    if isempty(fpath)
        warning('canlab_render_patterns:missing', ...
            'Image %s not found on path or in folder, skipping.', relname);
        continue
    end

    fprintf('\n=== %s ===\n  %s\n', label, fpath);

    try
        obj = fmri_data(fpath, 'noverbose');
    catch ME
        warning('canlab_render_patterns:fmridata', ...
            'Could not load %s as fmri_data: %s', fpath, ME.message);
        continue
    end

    if ~isempty(opts.thresh)
        try, obj = threshold(obj, opts.thresh); end %#ok<TRYNC>
    end

    if ismember('surface', opts.kinds)
        try
            fig = create_figure([label '_surface']);
            axis off;
            try
                surface(obj, 'foursurfaces_hcp', 'nolegend', 'nofigure');
            catch
                clf(fig); surface(obj, 'nolegend', 'nofigure');
            end
            drawnow;
            saveas(fig, fullfile(opts.savedir, [label '_surface.png']));
        catch ME
            warning('canlab_render_patterns:surface', ...
                '%s: surface failed (%s)', label, ME.message);
        end
    end

    if ismember('montage', opts.kinds)
        try
            fig = create_figure([label '_montage']);
            axis off;
            o2 = montage(obj, 'nolegend');
            try, title_montage(o2, 5, label); end %#ok<TRYNC>
            drawnow;
            saveas(fig, fullfile(opts.savedir, [label '_montage.png']));
        catch ME
            warning('canlab_render_patterns:montage', ...
                '%s: montage failed (%s)', label, ME.message);
        end
    end

    if ismember('isosurface', opts.kinds)
        try
            fig = figure('Name', [label '_isosurface']);
            han = isosurface(obj);
            if ~isempty(han)
                set(han, 'FaceAlpha', 0.5);
            end
            view(135, 20);
            try, lightFollowView; end %#ok<TRYNC>
            try, lightRestoreSingle; end %#ok<TRYNC>
            axis off; drawnow;
            saveas(fig, fullfile(opts.savedir, [label '_isosurface.png']));
        catch ME
            warning('canlab_render_patterns:isosurface', ...
                '%s: isosurface failed (%s)', label, ME.message);
        end
    end

end

fprintf('\nDone. PNGs written to: %s\n', opts.savedir);

end % function
