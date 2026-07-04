function run_all_multivariate(category_path)
% run_all_multivariate  Run visualize_contents.m for every study folder
% under a top-level category in Neuroimaging_Pattern_Masks.
%
% Usage:
%   run_all_multivariate                                              % default: Multivariate_signature_patterns
%   run_all_multivariate('CANlab_Meta_analysis_maps')                 % a different category
%   run_all_multivariate('/abs/path/to/Multivariate_signature_patterns')
%
% Re-runs all visualize_contents.m scripts to regenerate png_images/.
% Requires CanlabCore + SPM12 on path.

if nargin < 1 || isempty(category_path)
    category_path = 'Multivariate_signature_patterns';
end

% Resolve repo root from this file's location (docs/run_all_multivariate.m).
docs_dir = fileparts(mfilename('fullpath'));
if isempty(docs_dir), docs_dir = pwd; end
repo_root = fileparts(docs_dir);

if ~isfolder(category_path)
    candidate = fullfile(repo_root, category_path);
    if isfolder(candidate)
        category_path = candidate;
    else
        error('Could not find category folder: %s', category_path);
    end
end

% Make the shared helpers reachable from any folder's visualize_contents.m.
addpath(docs_dir);

set(0, 'DefaultFigureVisible', 'off');

d = dir(category_path);
folders = {};
for i = 1:numel(d)
    if ~d(i).isdir, continue, end
    if d(i).name(1) == '.', continue, end
    p = fullfile(category_path, d(i).name);
    if exist(fullfile(p, 'visualize_contents.m'), 'file')
        folders{end+1} = p; %#ok<AGROW>
    end
end

fprintf('Found %d folders with visualize_contents.m in %s\n', numel(folders), category_path);

for i = 1:numel(folders)
    fprintf('\n=== [%d/%d] %s ===\n', i, numel(folders), folders{i});
    try
        old_path = path;
        addpath(folders{i});
        cd(folders{i});
        visualize_contents();
        path(old_path);
    catch ME
        fprintf(2, '*** FAILED: %s\n%s\n', folders{i}, getReport(ME, 'extended'));
    end
    close all force;
end

fprintf('\nAll done.\n');

end
