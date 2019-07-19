% little script to make S1-back ROI as described in Lee et al 2018

% load any img in MNI space -- doesn't matter which
img = which('gray_matter_mask.img')
dat = fmri_data(img);
 
% peak coordinates given in paper (bilateral)
mm = [-18 -38 72; 18 -38 72];

% make 4mm radius sphere
indx = iimg_xyz2spheres(mm2voxel(mm,dat.volInfo),dat.volInfo.xyzlist,4);

% save back into img, view to confirm, and write it to file
dat.dat = indx;
orthviews(dat)
write(dat, 'fname', 'S1back_Lee_2018coords.nii')