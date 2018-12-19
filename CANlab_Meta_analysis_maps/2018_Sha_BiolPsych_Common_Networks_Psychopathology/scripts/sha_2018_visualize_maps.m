% Run this from the main directory above the "scripts" folder

% Display helper functions
% --------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

% Get file names
% --------------------------------------------------------

addpath('scripts')

mapdir = fullfile(pwd, 'MKDA_Maps');

imgs = dir(fullfile(mapdir, 'VBM', '*nii'));

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'DMN', 'DMN_Dec', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'DMN', 'DMN_Inc', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'DMN', 'DMN_Pool', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'FPN', 'FPN_Dec', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'FPN', 'FPN_Inc', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'FPN', 'FPN_Pool', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'SN', 'SN_Dec', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'SN', 'SN_Inc', '*nii'))];

imgs = [imgs; dir(fullfile(mapdir, 'SB_FC', 'SN', 'SN_Pool', '*nii'))];

% Load and display
% --------------------------------------------------------

for i = 1:length(imgs)
    
    img_to_load = fullfile(imgs(i).folder, imgs(i).name);
    
    printhdr(imgs(i).name)
    
    obj = fmri_data(img_to_load, 'noverbose');
    
    create_figure('montage'); axis off
    montage(obj);
    
    surface(obj, 'cutaway', 'ycut_mm', -35, 'noverbose')
    
    % montage(obj, 'full');
    %surface(obj, 'cutaway');
    
    drawnow, snapnow
    
end

