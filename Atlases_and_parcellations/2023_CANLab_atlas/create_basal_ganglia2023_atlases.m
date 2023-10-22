for scale = {'s1', 's4'}
    scale = scale{1};
    for space = {'fmriprep20','fsl6'}
        space = space{1};
        atlas_name = sprintf('bg_%s_%s',scale, space);
        atlas_obj = load_atlas(sprintf('tian_3t_%s_%s',scale,space));
        
        % Threshold some to clean up and avoid bleed-over
        % atlas_obj = threshold(atlas_obj, .6);
        
        atlas_obj = check_properties(atlas_obj);
        
        % I've asked the authors to give us probability maps around 10/13/23, but
        % until we get them we should assign an arbitrary probability to these all
        % so that they get along better with other atlases when merged.
        atlas_obj.probability_maps(atlas_obj.probability_maps == 1) = 0.8;
        
        atlas_obj = atlas_obj.select_atlas_subset(find(contains(atlas_obj.labels, {'PUT','CAU','NAc','GP'})));
        
        
        %% Check display
        % -----------------------------------------------------------------------
        
        % Display with unique colors for each region:
        orthviews(atlas_obj, 'unique');
         
        %% Convert to regions
        % -----------------------------------------------------------------------
        
         r = atlas2region(atlas_obj);
        
        %% Enforce some var types and compress
        
        atlas_obj = check_properties(atlas_obj);
        atlas_obj = remove_empty(atlas_obj);
        
        
        %% Save
        
        targetdir = what('2023_CANLab_atlas');
        cd(targetdir.path);
        
        savefile = sprintf('Basal_ganglia_%s_%s_combined_atlas_object.mat', scale, space);
        save(savefile, 'atlas_obj');
        
         %% save figure
        
        if dosave
           
            switch space
                case 'fmriprep20'
                    template = 'fmriprep20_template.nii.gz';
                case 'fsl6'
                    template = 'fsl6_hcp_template.nii.gz';
                otherwise
                    warning('Did not recognize space as one with an available underlay. Please update code accordingly if this is incorrect.');
            end
        
            o2 = canlab_results_fmridisplay([], 'full2', 'overlay', which(template));
            brighten(.6)
            
            o2 = montage(r, o2, 'wh_montages', 1:2);
            
            savedir = fullfile(pwd, 'png_images');
            if ~exist(savedir, 'dir'), mkdir(savedir); end
            
            scn_export_papersetup(600);
            savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
            saveas(gcf, savename);
        
            
        end
         
        %% Isosurface
        
        if dosave
            
            figure; han = isosurface(atlas_obj);
            
            cellfun(@(x1)set(x1,'FaceAlpha', .5), han)
            view(135, 20)
            lightFollowView;
            lightRestoreSingle
            axis off
            
            savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
            saveas(gcf, savename);
            
        end
    end
end