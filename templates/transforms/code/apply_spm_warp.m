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
    job = struct('comp', {{struct('def',{{[fname.folder, '/', fname.name]}}),...
                           struct('idbbvox', struct('vox', res(1:3)', ...
                                                    'bb', [-91, -126, -72; ...
                                                           90, 91, 100]))}},...
                 'out', {{struct('pull',struct('fnames', {mvg_img_vols},...
                                               'savedir',struct('savesrc',1),...
                                               'interp',interp, ...
                                               'mask', 1, ...
                                               'fwhm', [0,0,0],...
                                               'prefix','w'))}});

    spm_deformations(job);

    old_mvg_img = mvg_img;
    mvg_img_fname = dir(mvg_img);
    delete(old_mvg_img);
    wmvg_img = sprintf('%s/w%s', mvg_img_fname.folder, mvg_img_fname.name);

    fprintf('Warped data written to %s/w%s\n', mvg_img_fname.folder, mvg_img_fname.name);

    if ~isempty(post_affine_mat)
        M = csvread(post_affine_mat);
        wmvg_img_hdr = spm_vol(wmvg_img);
        for i = 1:length(wmvg_img_hdr)
            spm_get_space(sprintf('%s,%d',wmvg_img,i), M*wmvg_img_hdr(1).mat);
        end
        flags = struct('interp',interp, 'mask', 1, 'mean', 0, 'which', 2, 'wrap', zeros(3,1));
        spm_reslice({fxd_img, wmvg_img},flags);
        
        rwmvg_img = sprintf('%s/rw%s', mvg_img_fname.folder, mvg_img_fname.name);
        movefile(rwmvg_img, out_img);
        fprintf('Warped data moved to %s\n', out_img);
        delete(wmvg_img);
    else
        movefile(wmvg_img, out_img);
        fprintf('Warped data moved to %s\n', out_img);
    end

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
