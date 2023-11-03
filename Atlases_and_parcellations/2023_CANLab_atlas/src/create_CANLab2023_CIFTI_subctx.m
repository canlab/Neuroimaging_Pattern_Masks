function create_CANLab2023_CIFTI_subctx(SPACE,SCALE,res,atlas_obj)
    if res ~= 2
        error('We only have reference volumes for 2mm HCP spaces');
    end
    switch SPACE
        case 'MNI152NLin6Asym'
            cifti_mask = fmri_mask_image(which('hcp_cifti_subctx_labels.nii.gz'));
        case 'MNI152NLin2009cAsym'
            warning('Using projection from hcp grayordinate space into MNI152NLin2009cAsym. See sourcecode for details')
            % This code was developed using CIFTI templates derived from
            % HCP data which were spatially aligned into
            % MNI152NLin2009cAsym space. What we need are volumetric
            % templates taken directly from volumetric CIFTI segmentations 
            % of MNI152NLin2009cAsym aligned data.
            cifti_mask = fmri_mask_image(which('hcp_cifti_subctx_labels.nii.gz'));
        otherwise
            error('Unsupported space %s', SPACE)
    end

    this_dir = dir(which('create_CANLab2023_CIFTI_subctx.m'));
    %% output label files
    % we need these to match Glasser parcel labels to canlab atlas labels
    labels = atlas_obj.labels;
    fid = fopen(sprintf('%s/CANLab2023_%s_%s_%dmm_cifti_canlab_labels.txt', this_dir.folder, SCALE, SPACE, round(res)),'w+');
    for i = 1:num_regions(atlas_obj)
        fprintf(fid, '%d %s\n',i,labels{i});
    end
    fclose(fid);
    
    ctx_ind = contains(labels,'Ctx_');
    labels(ctx_ind) = cellfun(@(x1)regexprep(x1,'Ctx_(.*)_([LR])','$2_$1_ROI'),labels(ctx_ind),'UniformOutput',false);
    labels(ctx_ind) = cellfun(@(x1)(regexprep(regexprep(strrep(x1,'_','-'), '^([LR])-','$1_'),'-ROI$','_ROI')), labels(ctx_ind), 'UniformOutput', false);
    
    fid = fopen(sprintf('%s/CANLab2023_%s_%s_%dmm_cifti_glasser_labels.txt', this_dir.folder, SCALE, SPACE, round(res)),'w+');
    for i = 1:num_regions(atlas_obj)
        fprintf(fid, '%d %s\n',i,labels{i});
    end
    fclose(fid);
    
    n_labels = num_regions(atlas_obj);

    %% output volumetric file
    if any(atlas_obj.volInfo.mat ~= cifti_mask.volInfo.mat,'all')
        warning('Volumetric CIFTI mask is not aligned with atlas_obj, data will be resampled, which is not ideal.')
        % this data should all be in register by design, if not you need to
        % revisit your atlas generation script and create a version that's
        % in register. In particular small structures should be projected
        % and aligned with this space using the fewest number of
        % interpolation steps, and each resampling results in an
        % interpolation.

        atlas_obj = atlas_obj.resample_space(cifti_mask);
        if n_labels ~= num_regoins(atlas_obj)
            warning('Atlas parcels dropped during resampling to cifti mask, CIFTI labels ("keys") won''t be consistent with canlab atlas');
            delete(sprintf('%s/CANLab2023_%s_%s_%dmm_cifti_index_labels.txt', this_dir.folder, SCALE, SPACE, round(res)));
        end
    end

    atlas_obj = atlas_obj.threshold(0.2).apply_mask(cifti_mask);

    subctx_ind = find(~contains(atlas_obj.labels, 'Ctx'));
    subctx_lbls = atlas_obj.labels(subctx_ind);
    if any(contains(atlas_obj.labels, 'Ctx'))
        n_overlap = sum(atlas_obj.select_atlas_subset({'Ctx'}).remove_empty.dat ~= 0);
        warning('%d cortical voxels overlap with cifti_mask. Will remove. Check results for gaps',n_overlap)
        ind = find(~contains(atlas_obj.labels, 'Ctx'));
        atlas_obj = atlas_obj.select_atlas_subset(ind);
    end
    if any(~ismember(subctx_lbls, atlas_obj.labels))
        excluded = subctx_lbls(~ismember(subctx_lbls, atlas_obj.labels));
        warning('%s is presumed non-cortical but excluded from CIFTI mask.\n',excluded{:})
    end

    roi_vx_n = zeros(num_regions(atlas_obj),1);
    for i = 1:num_regions(atlas_obj)
        roi_vx_n(i) = sum(atlas_obj.dat == i);
    end

    if any(roi_vx_n < 10)
        warning('%s, retained for CIFTI volume, has few voxels (n<10)\n',atlas_obj.labels{roi_vx_n < 10});
    end
    atlas_obj.probability_maps = [];

    if round(res) ~= res
        warning('Rounding res label in filenmae to nearest integer. This may overwrite any existing labels of different resolution. Update this code to accomodate fractional resolutions more gracefully.');
    end
    atlas_obj.fullpath = sprintf('%s/CANLab2023_%s_%s_%dmm_cifti_vols.nii', this_dir.folder, SCALE, SPACE, round(res));

    atlas_obj.write('overwrite');
    gzip(atlas_obj.fullpath);
    delete(atlas_obj.fullpath);

    n_roi = length(atlas_obj.labels);
    cmap = round(255*colormap('lines'));
    cmap = repmat(cmap,ceil(n_roi/length(cmap)),1);

    fid = fopen(sprintf('%s/CANLab2023_%s_%s_%dmm_cifti_vols.txt', this_dir.folder, SCALE, SPACE, round(res)),'w+');
    for i = 1:n_roi
        fprintf(fid, [atlas_obj.labels{i}, '\n']);
        fprintf(fid, [int2str(i), ' ' num2str(cmap(i,1)), ' ', num2str(cmap(i,2)), ' ', num2str(cmap(i,3)), ' 255\n']);
    end
    fclose(fid);

    fprintf('To complete creation of CIFTI atlas please configure and run Atlases_and_parcellations/2023_CANLab_atlas/src/create_CANLab2023_atlas_cifti.sh');
end