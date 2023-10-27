function download_warpfield(from, to, format)
    % from and to should be formated as MNI152NLin6Asym or something standard like that. Options atm are
    % MNI152NLin6Asym
    % MNI152NLin2009cAsym
    % format should be ants, spm or fsl

    this_dir = dir(which('download_warpfield.m'));

    transform = [from, '-', to];
    switch transform
        % in the code below we download to a tempname in case the download
        % is interrupted. We don't want file existence checks to flag a partial
        % download as a valid file.
        case 'MNI152NLin6Asym-MNI152NLin2009cAsym'
            %hosted on bogpetre@gmail.com's figshare
            file = tempname;
            websave(file, 'https://figshare.com/ndownloader/files/42771256');
            movefile(file, [this_dir.folder, '/', format, '/y_01_fsl_to_fmriprep_DisplacementFieldTransform.nii']);

            file = tempname;
            websave(file, 'https://figshare.com/ndownloader/files/42771256');
            movefile(file, [this_dir.folder, '/', format, '/y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii']);
        case 'MNI152NLin2009cAsym-MNI152NLin6Asym'
            file = tempname;
            websave(file, 'https://figshare.com/ndownloader/files/42771259');
            movefile(file, [this_dir. folder, '/', format, '/y_00_fmriprep_to_fsl_DisplacementFieldTransform.nii']);

            file = tempname;
            websave(file, 'https://figshare.com/ndownloader/files/42771256');
            movefile(file, [this_dir.folder, '/', format, '/y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii']);
        otherwise
            error('Could not find transform from %s to %s in %s format', from, to, format);
    end
end
