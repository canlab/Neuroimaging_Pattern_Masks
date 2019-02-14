mydir = what('2018_Sha_BiolPsych_Common_Networks_Psychopathology');

if isempty(mydir.path), error('Cannot find main directory on path.'), end

cd(mydir.path)

scriptsdir = fullfile(mydir.path, 'scripts');
addpath(scriptsdir)

pubfilename = fullfile(scriptsdir, 'sha_2018_visualize_maps');

pubdir = fullfile(mydir.path, 'html_summary_report');
if ~exist(pubdir, 'dir'), mkdir(pubdir), end

% ------------------------------------------------------------------------
%pubfilename = fullfile(pubdir,

p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 1600, ...
    'format', 'html', 'outputDir', pubdir, ...
    'showCode', false);

publish(pubfilename, p)

close all
