function visualize_contents(varargin)
% visualize_contents  Render Neurosynth 100-topic FDR maps into png_images/.
% See contents_description.md.
%
% Optional name/value:
%   'whichtopics' : numeric IDs to render (e.g. [13 17 60 61 97]).
%                   Default [] = render all 50 topics (slow!).
%   'mode'        : 'pAgF' | 'pFgA' | 'both' (default 'pAgF' to keep
%                   run-time reasonable).

p = inputParser;
p.addParameter('whichtopics', [], @isnumeric);
p.addParameter('mode', 'pAgF', @(s) any(strcmpi(s, {'pAgF','pFgA','both'})));
p.parse(varargin{:});

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

modes = {};
if any(strcmpi(p.Results.mode, {'pAgF','both'})), modes{end+1} = 'pAgF'; end
if any(strcmpi(p.Results.mode, {'pFgA','both'})), modes{end+1} = 'pFgA'; end

imgs = {};
for m = 1:numel(modes)
    pat = sprintf('v4-topics-100_*_%s_z_FDR_0.01.nii.gz', modes{m});
    ff = dir(fullfile(this_dir, pat));
    for i = 1:numel(ff)
        % Topic ID is the first underscore-delimited number after "v4-topics-100_".
        tok = regexp(ff(i).name, '^v4-topics-100_(\d+)_', 'tokens', 'once');
        if isempty(tok), continue, end
        tid = str2double(tok{1});
        if ~isempty(p.Results.whichtopics) && ~ismember(tid, p.Results.whichtopics)
            continue
        end
        label = regexprep(ff(i).name, '\.nii\.gz$', '');
        label = regexprep(label, '_z_FDR_0\.01$', '');
        imgs(end+1, :) = {label, ff(i).name}; %#ok<AGROW>
    end
end

if isempty(imgs)
    warning('visualize_contents:notopics', 'No matching topic files found.');
    return
end

canlab_render_patterns(imgs);

end
