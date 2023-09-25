function dilated_atlas = dilate(img, mask)
    % takes an atlas object img and dilates it by searching for the nearest
    % label for each voxel in mask. Runs much faster if you resample img to
    % mask first, but it will be less accurate

    img = img.replace_empty();

    % dilation does not cross midline    
    % split into halves
    xyz = img.volInfo.mat*[img.volInfo.xyzlist'; ones(1,size(img.dat,1))];
    imgL = img;
    imgL.dat(xyz(1,:) >= 0) = 0;
    imgR = img;
    imgR.dat(xyz(1,:) < 0) = 0;

    img = img.remove_empty();
    imgL = imgL.remove_empty();
    imgR = imgR.remove_empty();

    nearest = zeros(size(mask.volInfo.xyzlist,1),1);
    ind = find(mask.dat);
    target_mm = mask.volInfo.mat*[mask.volInfo.xyzlist(ind,:)';ones(1,length(ind))];
    target_mm = target_mm(1:3,:);
    for i = 1:length(ind)
        if target_mm(1,i) < 0
            this_img = imgL;
        else
            this_img = imgR;
        end
        this_ind = ind(i);
        this_region = this_img.find_closest_region(target_mm(:,i)).region_number;
        nearest(this_ind) = this_region;
    end
    
    dilated_atlas = img.resample_space(mask).replace_empty();

    dilated_atlas.dat = nearest;
    dilated_atlas.probability_maps(~mask.dat,:) = 0;
end