fname = which('MNI152NLin2009cAsym_margulies_grad1.nii.gz')
obj = fmri_data(fname);

create_figure('montage'); axis off; montage(obj);
saveas(gcf, 'Margulies_2009asym_montage.png');

create_figure('surface'); axis off; surface(obj, 'foursurfaces');
saveas(gcf, 'Margulies_2009asym_foursurfaces.png');