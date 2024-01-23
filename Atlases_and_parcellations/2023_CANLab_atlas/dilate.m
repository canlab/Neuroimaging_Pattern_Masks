function dilated_atlas = dilate(img, mask)
    % takes an atlas object img and dilates it by searching for the nearest
    % label for each voxel in mask. Runs much faster if you resample img to
    % mask first, but it will be less accurate
    %
    % Note, this script will not dilate areas across the midline, but if
    % the input is bilateral it may already cross the midline, in which
    % case a region that is majority say right lateralized may have a voxel
    % or two on the midline or in the left hemisphere, in which case that
    % label may end up impact dilations in the contralateral hemisphere
    % regardless. If you need this script to produce lateralized regions
    % it's best to manually separate L and R side regions and dilate each
    % separately before recombining.

    mask = mask.resample_space(img,'nearest');

    img = img.replace_empty();

    % dilation does not cross midline    
    % split into halves, assigning midline to x+ side
    xyz = img.volInfo.mat*[img.volInfo.xyzlist'; ones(1,size(img.dat,1))];
    imgL = img;
    imgL.dat(xyz(1,:) >= 0) = 0;
    imgR = img;
    imgR.dat(xyz(1,:) < 0) = 0;

    imgL = imgL.remove_empty();
    imgR = imgR.remove_empty();

    nearest = zeros(size(mask.volInfo.xyzlist,1),1);
    new_voxels = find(mask.apply_mask(fmri_mask_image(img),'invert').replace_empty.dat); % only consider areas not already in img
    target_mm = mask.volInfo.mat*[mask.volInfo.xyzlist(new_voxels,:)';ones(1,length(new_voxels))];
    target_mm = target_mm(1:3,:);
    region_found = zeros(length(new_voxels),1);
    parfor i = 1:length(new_voxels)
        if target_mm(1,i) < 0
            this_img = imgL;
        else
            this_img = imgR;
        end
        region_found(i) = this_img.find_closest_region(target_mm(:,i)).region_number;
    end

    for i = 1:length(new_voxels)
        this_ind = new_voxels(i);
        nearest(this_ind) = region_found(i);
    end
    
    dilated_atlas = img.replace_empty();

    dilated_atlas.dat(new_voxels) = nearest(new_voxels);
    if ~isempty(dilated_atlas.probability_maps)
        dilated_atlas.probability_maps(~mask.dat,:) = 0;
        for i = 1:length(new_voxels)
            roi = dilated_atlas.dat(new_voxels(i));
            % the value we assign here is arbitrary, but we know this
            % expansion is artificial so let's cede the region to whatever
            % competition we may encounter from other ROIs
            dilated_atlas.probability_maps(new_voxels(i),roi) = 0.2;
        end
    end
    
    % drop any null ROIs
    new_voxels = zeros(1,length(dilated_atlas.labels));
    new_voxels(unique(dilated_atlas.dat(dilated_atlas.dat ~=0))) = true;

    remove = ~new_voxels & ~any(dilated_atlas.probability_maps);
    keep = find(~remove);
    remove = find(remove);
    if any(remove)
        for r = 1:length(remove)
            warning('Dropping ROI %d',remove(r))
        end
    end

    dilated_atlas = dilated_atlas.select_atlas_subset(keep);
end