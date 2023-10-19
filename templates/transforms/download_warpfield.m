function download_warpfield(from, to, format)
    % from and to should be formated as MNI152NLin6Asym or something standard like that. Options atm are
    % MNI152NLin6Asym
    % MNI152NLin2009cAsym
    $ format should be ants, spm or fsl

    this_dir = dir(which('download_spm_transforms.m'));

    transform = [from, '-', to];
    switch :
        case 'MNI152NLin6Asym-MNI152NLin2009cAsym':
            # hosted on bogpetre@gmail.com's figshare
            websave([this_dir.folder, '/', format, '/y_01_fsl_to_fmriprep_DisplacementFieldTransform.nii'], ...
                'https://figshare.com/ndownloader/files/42771256');
        case 'MNI152NLin2009cAsym-MNI152NLin6Asym':
            websave([this_dir. folder, '/', format, '/y_00_fmriprep_to_fsl_DisplacementFieldTransform.nii'], ...
                'https://figshare.com/ndownloader/files/42771259');
        otherwise
            error('Could not find transform from %s to %s in %s format', from, to, format);
    end
end
