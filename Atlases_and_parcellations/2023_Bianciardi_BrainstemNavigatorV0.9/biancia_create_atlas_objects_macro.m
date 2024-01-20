% For most atlases in this repo the atlas generation script serves as
% documentation of how an atlas was created. In the case of the Bianciardi
% Brainstem Navigator atlas we can't distribute the atlas itself, so
% instead of a record of how an atlas was created we have scripts that
% download and create the atlas dynamically for each person who invokes the
% load_atlas() command for one of these atlases or runs this script. Unlike
% the load_atlas path, this script will generate all the usual accompanying
% files that most atlases have, like the png_images and the regions files.
% It also serves as a useful oneshot script to download and assemble all
% variations of atlas in one shot, rather than waiting for the first time
% you invoke load_atlas() on them.

clear all; close all;

addpath('~/.matlab/spm/spm12');

addpath(genpath('~/.matlab/canlab/CanlabCore'))
addpath(genpath('~/.matlab/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('~/.matlab/canlab/MasksPrivate'))

for space = {'MNI152NLin6Asym', 'MNI152NLin6Asym_2mm', 'MNI152NLin2009cAsym', 'MNI152NLin2009cAsym_2mm'}
    space = space{1};
    switch space
        case {'MNI152NLin6Asym','MNI152NLin6Asym_2mm'}
            alias = 'fsl6';
        case {'MNI152NLin2009cAsym','MNI152NLin2009cAsym_2mm'}
            alias = 'fmriprep20';
        otherwise
            alias = space;
    end

    bianciaAtlas = bianciardi_create_atlas_obj(space);

    atlas_name = bianciaAtlas.atlas_name;
    % Run this from the directory containing the atlas files
    % -----------------------------------------------------------------------
    dosave = true;
    
    % Check display
    % -----------------------------------------------------------------------
    
    % Display with unique colors for each region:
    
    overlayargs = {};
    if ~isempty(dir(which([alias, '_template.nii.gz'])))
        overlayargs = {'overlay', which([alias, '_template.nii.gz'])};
    elseif ~isempty(dir(which([alias, '_template.nii'])))
        overlayargs = {'overlay', which([alias, '_template.nii.gz'])};
    end
    
    orthviews(bianciaAtlas, 'unique', overlayargs{:});
    figure;
    % Convert to regions
    % -----------------------------------------------------------------------
    
    r = atlas2region(bianciaAtlas);
    
    % Display on montage (colors may not be the same!):
    % montage(r);
     
    %% save figure
    if dosave
        o2 = canlab_results_fmridisplay([], 'full2', overlayargs{:});
        brighten(.6)
        
        o2 = montage(r, o2);
        
        savedir = fullfile(pwd, 'png_images');
        if ~exist(savedir, 'dir'), mkdir(savedir); end
        
        scn_export_papersetup(600);
        savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
        saveas(gcf, savename);
    
        
    end
    
    
    %% write - this writes only the label image
    
    if dosave
        
        savename = sprintf('%s_atlas_regions.img', atlas_name);
        bianciaAtlas.fullpath = fullfile(pwd, savename);
        write(bianciaAtlas,'overwrite');
        
    end
    
    %% Turn regions into separate list of names, for canlab_load_ROI
    % which loads regions by name from mat files.
    
    clear region_names
    
    for i = 1:length(r)
        
        eval([bianciaAtlas.labels{i} ' = r(i);']);
        
        region_names{i} = r(i).shorttitle;
        
    end
    
    savename = sprintf('%s_atlas_regions.mat', atlas_name);
    save([pwd, '/' savename], 'r', 'region_names', bianciaAtlas.labels{:});
    
    %%
    if dosave
        
        figure; han = isosurface(bianciaAtlas);
        
        arrayfun(@(x1)set(x1,'FaceAlpha', .5), han)
        view(135, 20)
        lightFollowView;
        lightRestoreSingle
        axis off
        
        pos = get(gcf,'position');
        set(gcf,'position',[pos(1:2), 570, 674]);
    
        savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
        saveas(gcf, savename);
        
    end

end