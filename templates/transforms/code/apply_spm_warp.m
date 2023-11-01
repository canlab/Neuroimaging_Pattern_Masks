% this script applies transforms using spm and custom code (I don't know how SPM applies affine transforms, it's not well
% documented, so I implemented something myself built off SPM utils).

function apply_spm_warp(mvg_img0, fxd_img0, pre_affine_mat, warp_img, post_affine_mat, out_img, interp)
    % all arguments are paths except for interp which is an integer corresponding to the following
    % 0 - 'nearest neighbour'
    % 1 - 'trilinear'
    % 2 - 'spline2' - second degree b-spline
    % 3 - 'spline3' - third degree b-spline
    % 4 - 'spline4' - fourth degree b-spline
    % 5 - 'spline5' - fifth degree b-spline
    % 6 - 'spline6' - sixth degree b-spline
    % 7 - 'spline7' - seventh degree b-spline
    %
    % this function is designed to work with the prepared transformations, not arbitrary ones.
    % I couldn't find helpful documentation for SPM and can't make any guarantees that this will work with any other
    % transforms, but if you do want to adapt it have a look at the bb (bounding box) option below. This should be
    % asigned dynamically based on the input and/or output image at the very least.
    %
    % This script has only been tested for alignment of 3D templates to one another.

    % make a copy so we don't modify original if we pretransform or
    % posttransform
    fname = dir(mvg_img0);
    mvg_img = sprintf('%s.nii', tempname);
    copyfile(mvg_img0, mvg_img);
    
    fname = dir(fxd_img0);
    fxd_img = sprintf('%s.nii', tempname);
    copyfile(fxd_img0, fxd_img);

    mvg_img_hdr = spm_vol(mvg_img);
    n_imgs = length(mvg_img_hdr);
    fxd_img_hdr = spm_vol(fxd_img0);


    if ~isempty(pre_affine_mat)
        M = csvread(pre_affine_mat);
        for i = 1:length(mvg_img_hdr)
            spm_get_space(sprintf('%s,%d',mvg_img,i), M*mvg_img_hdr(1).mat);
        end
        
        res = diag(fxd_img_hdr.mat);
    else
        res = diag(mvg_img_hdr(1).mat);
    end

    


    zippedWarp = false;
    if contains(warp_img,'nii.gz')
        zippedWarp = true;
        warp_img_new = [tempname, '.nii.gz'];
        copyfile(warp_img, warp_img_new)
        gunzip(warp_img_new)
        warp_img = warp_img_new;
    end
    fname = dir(warp_img);
    mvg_img_vols = cell(length(mvg_img_hdr),1);
    for i = 1:length(mvg_img_hdr)
        mvg_img_vols{i} = [mvg_img, ',', int2str(i)];
    end
    %{
    job = struct('comp', {{struct('def',{{[fname.folder, '/', fname.name]}}),...
                           struct('idbbvox', struct('vox', res(1:3)', ...
                                                    'bb', [-91, -126, -72; ...
                                                           90, 91, 100]))}},...
    %}
    bbox = spm_get_bbox(fxd_img_hdr);
    job = struct('comp', {{struct('def',{{[fname.folder, '/', fname.name]}}),...
                           struct('idbbvox', struct('vox', res(1:3)', ...
                                                    'bb', bbox))}},...
                 'out', {{struct('pull',struct('fnames', {mvg_img_vols},...
                                               'savedir',struct('savesrc',1),...
                                               'interp',interp, ...
                                               'mask', 0, ...
                                               'fwhm', [0,0,0],...
                                               'prefix','w'))}});

    spm_deformations(job);
    mvg_img_fname = dir(mvg_img);
    wmvg_img = sprintf('%s/w%s', mvg_img_fname.folder, mvg_img_fname.name);

    old_mvg_img = mvg_img;
    delete(old_mvg_img);

    fprintf('Warped data written to %s/w%s\n', mvg_img_fname.folder, mvg_img_fname.name);

    if ~isempty(post_affine_mat)
        M = csvread(post_affine_mat);
        final_img_hdr = spm_vol(wmvg_img);
        for i = 1:length(final_img_hdr)
            spm_get_space(sprintf('%s,%d',wmvg_img,i), M*final_img_hdr(1).mat);
        end
        flags = struct('interp',interp, 'mask', 0, 'mean', 0, 'which', 2, 'wrap', zeros(3,1));
        spm_reslice({fxd_img, wmvg_img},flags);
        
        final_img = sprintf('%s/rw%s', mvg_img_fname.folder, mvg_img_fname.name);
    else
        final_img = wmvg_img;
    end

    % reorient to reference if needed
    final_img_hdr = spm_vol(final_img);
    for i = 1:length(final_img_hdr)
        assert(all(final_img_hdr(1).mat == final_img_hdr(i).mat,'all'),'Transformed images have different header matrices across volumes. This shouldn''t happen');
        assert(all(final_img_hdr(1).dim == final_img_hdr(i).dim,'all'), 'Transformed images have different volume sizes. This shouldn''t happen');
    end
    
    reorient_mat = final_img_hdr(1).mat\fxd_img_hdr(1).mat;
    if any(reorient_mat ~= eye(4),'all')
        assert(all(abs(diag(reorient_mat)) == ones(4,1)),'Image was not correctly resampled');

        % flip the data and header matrix
        V = spm_read_vols(final_img_hdr);
        [X,Y,Z] = ndgrid(1:size(V,1), 1:size(V,2), 1:size(V,3));
        XYZ = [X(:),Y(:),Z(:)]';
        XYZ_new = reorient_mat*[XYZ; ones(1,numel(X))];

        for i = 1:length(final_img_hdr)
            V0 = V(:,:,:,i);
            V1 = reshape(V0(sub2ind(size(V0), XYZ_new(1,:), XYZ_new(2,:), XYZ_new(3,:))), size(V0));

            final_img_hdr(i).mat = fxd_img_hdr(1).mat;
            spm_write_vol(final_img_hdr(i),V1);
        end

        % write it out
    end

    movefile(final_img, out_img);
    fprintf('Warped data moved to %s\n', out_img);


    if ~isempty(dir(mvg_img))
        delete(mvg_img)
    end

    if ~isempty(dir(fxd_img))
        delete(fxd_img)
    end

    if zippedWarp
        delete(warp_img)
    end
end
