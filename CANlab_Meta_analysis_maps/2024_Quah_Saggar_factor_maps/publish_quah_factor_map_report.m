pubdir = fullfile(pwd, 'published_output');
if ~exist(pubdir, 'dir'), mkdir(pubdir), end

pubfilename = ['quah_factor_maps' scn_get_datetime];

p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 1600, ...
    'format', 'html', 'outputDir', fullfile(pubdir, pubfilename), 'showCode', false);

% ------------------------------------------------------------------------
% Run and report status
% ------------------------------------------------------------------------

publish('quah_factor_map_montages.m', p)

