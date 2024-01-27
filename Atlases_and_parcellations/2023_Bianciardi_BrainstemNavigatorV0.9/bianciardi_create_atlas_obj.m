function bianciaAtlas = bianciardi_create_atlas_obj(space)
    % build the bianciardi atlas in the specified space either the coarse
    % parcellation or the fine parcellation
    % space - 'MNI152NLin6Asym', 'MNI152NLin6Asym_2mm', 'MNI152NLin2009cAsym', or 'MNI152NLin2009cAsym_2mm'
    %
    % Note that this atlas has many tiny regions. It really helps to
    % minimize the number of interpolation steps, hence the 2mm versions of
    % atlas included here. It's much better to incorporate transformations
    % directly into 2mm space rather than transform into 1mm and then
    % resample using canlab tools. Take note if using any other spaces
    % besides those already included.
    this_dir = dir(which('bianciardi_create_atlas_obj.m'));

    switch space
        case 'MNI152NLin6Asym'
            space_lbl = 'fsl6';
            ref_file_base = 'MNI152NLin6Asym_T1_1mm';
        case 'MNI152NLin6Asym_2mm'
            space_lbl = 'fsl6_2mm';
            ref_file_base = 'MNI152NLin6Asym_T1_1mm';
        case 'MNI152NLin2009cAsym'
            space_lbl = 'fmriprep20';
            ref_file_base = 'MNI152NLin2009cAsym_T1_1mm';
        case 'MNI152NLin2009cAsym_2mm'
            space_lbl = 'fmriprep20_2mm';
            ref_file_base = 'MNI152NLin2009cAsym_T1_2mm';
        otherwise
            error('Unrecognized space %s',space);
    end
    SCALE='fine';
    ref_file = which([ref_file_base '.nii.gz']);
    if isempty(dir(ref_file)), ref_file = which([ref_file_base '.nii']); end
    if isempty(dir(ref_file)), error('Could not locate reference file for %s', space); end

    atlas_name = sprintf('bianciardi_%s', space);
    space_description = space;
    references = char({'García-Gomar MG, Videnovic A, Singh K, Stauder M, Lewis LD, Wald LL, Rosen BR, Bianciardi M. Disruption of brainstem structural connectivity in RBD using 7 Tesla MRI. Mov Disord. 2021 Dec 29. doi: 10.1002/mds.28895. Online ahead of print. PMID: 34964520',...
        'Singh K, García-Gomar MG, Bianciardi M. Probabilistic Atlas of the Mesencephalic Reticular Formation, Isthmic Reticular Formation, Microcellular Tegmental Nucleus, Ventral Tegmental Area Nucleus Complex, and Caudal-Rostral Linear Raphe Nucleus Complex in Living Humans from 7 Tesla Magnetic Resonance Imaging. Brain Connect. 2021 Oct;11(8):613-623. doi: 10.1089/brain.2020.0975. Epub 2021 biancian 17. PMID: 33926237.',...
        'Singh K, Indovina I, Augustinack JC, Nestor K, García-Gomar MG, Staab JP, Bianciardi M. Probabilistic Template of the Lateral Parabrachial Nucleus, Medial Parabrachial Nucleus, Vestibular Nuclei Complex, and Medullary Viscero-Sensory-Motor Nuclei Complex in Living Humans From 7 Tesla MRI. Front Neurosci. 2020 Jan 23;13:1425. doi: 10.3389/fnins.2019.01425. PMID: 32038134; PMCID: PMC6989551.',...
        'García-Gomar MG, Strong C, Toschi N, Singh K, Rosen BR, Wald LL, Bianciardi M. In vivo Probabilistic Structural Atlas of the Inferior and Superior Colliculi, Medial and Lateral Geniculate Nuclei and Superior Olivary Complex in Humans Based on 7 Tesla MRI. Front Neurosci. 2019 Aug 7;13:764. doi: 10.3389/fnins.2019.00764. PMID: 31440122; PMCID: PMC6694208.',...
        'Bianciardi M, Strong C, Toschi N, Edlow BL, Fischl B, Brown EN, Rosen BR, Wald LL. A probabilistic template of human mesopontine tegmental nuclei from in vivo 7T MRI. Neuroimage. 2018 Apr 15;170:222-230. doi: 10.1016/j.neuroimage.2017.04.070. Epub 2017 May 3. PMID: 28476663; PMCID: PMC5670016.',...
        'Bianciardi M, Toschi N, Edlow BL, Eichner C, Setsompop K, Polimeni JR, Brown EN, Kinney HC, Rosen BR, Wald LL. Toward an In Vivo Neuroimaging Template of Human Brainstem Nuclei of the Ascending Arousal, Autonomic, and Motor Systems. Brain Connect. 2015 Dec;5(10):597-607. doi: 10.1089/brain.2015.0347. Epub 2015 Aug 11. PMID: 26066023; PMCID: PMC4684653.'});
    
    % imort atlas file in MNI152NLin2009cAsym space
    % we bianciast use this as a stand in template that we'll modify later, since
    % this is the wrong space
    bianciaTbl = readtable(which('bianciardi_fine_labels.csv'));
    bianciaTbl_coarse = readtable(which('bianciardi_coarse_labels.csv'));
    labels = {};
    labels_2 = {};
    labels_3 = {};
    labels_4 = {};
    label_descriptions = {};
    for i = 1:height(bianciaTbl)
        labels{end+1} = bianciaTbl.left{i};
        label_descriptions{end+1} = bianciaTbl.full_label{i};
        labels_2{end+1} = bianciaTbl_coarse.left{i};
        labels_3{end+1} = bianciaTbl.Structure{i};
        labels_4{end+1} = num2str(bianciaTbl.Ref(i));
        if ~isempty(bianciaTbl.right{i})
            label_descriptions{end} = ['Left ',label_descriptions{end}];
            labels{end+1} = bianciaTbl.right{i};
            label_descriptions{end+1} = ['Right ' bianciaTbl.full_label{i}];    
            labels_2{end+1} = bianciaTbl_coarse.right{i};
            labels_3{end+1} = bianciaTbl.Structure{i};
            labels_4{end+1} = num2str(bianciaTbl.Ref(i));
        end
    end
    % fix capitalization
    for i = 1:length(label_descriptions), label_descriptions{i}(2:end) = lower(label_descriptions{i}(2:end)); end
    
    areaFile = get_area_file(labels, labels_3, this_dir);

    %% download files
    if any(cellfun(@isempty,areaFile))
        % Download necessary files

        fprintf('Bianciardi atlas source files not found in Matlab path.\n');        

        x = input('Would you like to attempt to download source data from the atlas developers?? y/n\n','s');

        if ~ismember(x,{'y','Y','yes','Yes','YES'})
            error('Could not find bianciardi atlas source files. This atlas has a restrictive license and cannot be distributed directly by this repo. It must be downloaded from the the Brainstem Imaging Lab (https://brainstemimaginglab.martinos.org/research/).');
        end

        %% License prompt
        fid = fopen([this_dir.folder, '/Copyright.txt'],'r');
        disp(fscanf(fid,'%c'));
        fclose(fid);
        
        x = input('Do you agree to this license?? y/n\n','s');

        if ~ismember(x,{'y','Y','yes','Yes','YES'})
            error('You may not use this atlas via this interface.');
        end

        %% download data
        outfile = fullfile(this_dir.folder, 'brainstemnavigatorv09.zip');
        fprintf('Saving file to %s\n',outfile);
        websave(outfile, 'https://www.nitrc.org/frs/download.php/12427/BrainstemNavigator.zip/','i_agree','1','release_id','4544','download_now', '1');
        fprintf('Extracting file to %s\n',this_dir.folder)
        unzip(outfile,this_dir.folder);
        delete(outfile);
        rmdir([this_dir.folder, '/BrainstemNavigator/0.9/1a.BrainstemNucleiAtlas_IIT'],'s')
        rmdir([this_dir.folder, '/BrainstemNavigator/0.9/1b.DiencephalicNucleiAtlas_IIT'],'s')
        rmdir([this_dir.folder, '/BrainstemNavigator/0.9/1c.Templates_IIT'],'s')

        areaFile = get_area_file(labels, labels_2, this_dir);
    end

    % import recently downloaded files
    % this version is faster but requires a lot of temp disk space
    %{
    area_path = cell(length(areaFile),1);
    for i=1:length(areaFile)
        this_file = areaFile{i};
        %area_path{i} = sprintf('%s/%s',this_file.folder, this_file.name);
    end
    original_parcel = fmri_data(area_path);
    %}
    area_path = cell(length(areaFile),1);
    for i=1:length(areaFile)
        this_file = areaFile{i};
        area_data{i} = fmri_data(sprintf('%s/%s',this_file.folder, this_file.name), 'noverbose');
    end
    original_parcel = cat(area_data{:});

    %% transform files
    switch space
        % identify necessary transformation files from atlas native
        % space (MNI152NLin6Asym) to desired space
        case {'MNI152NLin2009cAsym', 'MNI152NLin2009cAsym_2mm'}
            premat = which('00_fsl_to_fmriprep_subctx_AffineTransform.csv');
            warp = which('y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii');
            postmat = [];
            if isempty(warp)
                % this should be updated in the future to download the
                % transformation matrices from figshare
                fprintf('Could not find transformation matrices to project atlas into desired space. Downloading them from figshare.');
                
                download_warpfield('MNI152NLin6Asym', 'MNI152NLin2009cAsym', 'spm');
                
                warp = which('y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii');
            end
        case {'MNI152NLin6Asym','MNI152NLin6Asym_2mm'}
            premat = [];
            postmat = [];
            warp = [];
        otherwise
            error('Could not find transformation matrices to project atlas into desired space');
    end

    switch space
        case {'MNI152NLin2009cAsym','MNI152NLin2009cAsym_2mm'}
            % project data to necessary space if needed
            fprintf('Transforming from MNI152NLin6Asym to %s space\n', space);
    
            original_parcel.fullpath = [tempname, '.nii'];
            original_parcel.write('overwrite');
    
            aligned = [tempname, '.nii'];
            % align using precomputed warp fields and trilinear
            % interpolation
            apply_spm_warp(original_parcel.fullpath, ref_file,...
                premat,warp,postmat,...
                aligned, 1);
            delete(original_parcel.fullpath);
            pdata = fmri_data(aligned);
            delete(aligned)
        case 'MNI152NLin6Asym'
            % do nothing
            pdata = original_parcel;
        case 'MNI152NLin6Asym_2mm'
            pdata = original_parcel.resample_space(ref_file);
        otherwise
            error('Unrecognized space %s',space);        
    end

    % super low value errors are likely resampling issues. Not a problem
    % for FSL6 space, but a problem for any other spaces we might resample
    % to in other versions of this script. The issue with renormalization
    % is that there are boundary regions and we want to allow for a voxel
    % to be either one region or no region. i.e. if we have a p(A) = 0.3
    % and p(~A) = 0 we don't want to convert that into p(A) = 1 after
    % renormalization. We want to recognize that p(null) = 0.7. 
    pdata.dat(pdata.dat < 1e-4) = 0;

    %% renormalize probabilities
    % the probabilities don't add to one. There are two causes. First
    % different studies produced different sets of these probability maps.
    % Most studies have probabilities that sum to 1, but across studies
    % there are boundary overlaps. The second cause is specific to certain
    % regions in certain studies that overlap for some inexplicable reasons
    % despite the fact that they shouldn't. These are as follows, indexing
    % studies by publication year,
    % 2015 - 61 DR and PAG voxels overlap
    % 2021 - 1 voxels of L_VTA_PBP overlaps with L_mRta
    % In both cases we can just renormalize to fix the problem. 
    total_p = sum(pdata.dat,2);
    pdata.dat(total_p~=0,:) = pdata.dat(total_p~=0,:)./total_p(total_p~=0);
    
    %% build atlas
    
    % combine data with labels
    bianciaAtlas = atlas(pdata, ...
        'atlas_name', atlas_name,...
        'labels',labels, ...
        'label_descriptions', label_descriptions(:), ... 
        'labels_2', labels_2, ...
        'labels_3', labels_3, ...
        'labels_4', labels_4, ...
        'space_description', space_description, ...
        'references',references, 'noverbose');
    bianciaAtlas = bianciaAtlas.replace_empty();

    % prefix laterality for consistency with other atlases
    for i = 1:length(labels), labels{i} = regexprep(labels{i},'(.*)_l$','L_$1'); end
    for i = 1:length(labels), labels{i} = regexprep(labels{i},'(.*)_r$','R_$1'); end

    bianciaAtlas.labels = labels;

    % add b labels to raphe
    raphe_ind = contains(bianciaAtlas.label_descriptions,'pallidus');
    bianciaAtlas.labels{raphe_ind}(end+1:end+3) = '_B1';
    raphe_ind = contains(bianciaAtlas.label_descriptions,'obscurus');
    bianciaAtlas.labels{raphe_ind}(end+1:end+3) = '_B2';
    raphe_ind = contains(bianciaAtlas.label_descriptions,'magnus');
    bianciaAtlas.labels{raphe_ind}(end+1:end+3) = '_B3';
    raphe_ind = find(contains(bianciaAtlas.labels,'L_PnO_PnC'));
    bianciaAtlas.labels{raphe_ind}(end+1:end+3) = '_B5';
    raphe_ind = find(contains(bianciaAtlas.labels,'R_PnO_PnC'));
    bianciaAtlas.labels{raphe_ind}(end+1:end+3) = '_B5';
    raphe_ind = contains(bianciaAtlas.label_descriptions,'Dorsal');
    bianciaAtlas.labels{raphe_ind}(end+1:end+3) = '_B7'; 
    raphe_ind = find(contains(bianciaAtlas.labels,'Mn'));
    for ind = raphe_ind
        bianciaAtlas.labels{ind}(end+1:end+6) = '_B6_B8'; 
    end

    bianciaAtlas = threshold(bianciaAtlas, .01);
    bianciaAtlas.probability_maps = sparse(double(bianciaAtlas.probability_maps));
    timestamp = posixtime(datetime('Now'));
    bianciaAtlas.additional_info = struct('creation_date', {posixtime(datetime('Now'))});

    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save([this_dir.folder, '/' savename], 'bianciaAtlas');
    
    % we can't upload the mat file to github due to licensing issues, but
    % we can upload a timestamp that will flag out of date versions and
    % cause other uesrs to recreate the atlas when appropriate.
    fid = fopen(sprintf('%s/%s_atlas_object.latest', this_dir.folder, bianciaAtlas.atlas_name), 'w+');
    fprintf(fid,'%f',timestamp);
    fclose(fid);
end

function areaFile = get_area_file(labels, labels_2, parentDir)
    [areaFile, areaName] = deal(cell(length(labels),1));
    for i = 1:length(labels)
        areaName{i} = labels{i};
        switch areaName{i}
            case {'isRt_l', 'isRt_r','mRt_l'}
                areaName{i} = strrep(areaName{i},'t','T');
            case {'mRta_l','mRta_r'}
                areaName{i} = strrep(areaName{i},'t','T');
                areaName{i} = strrep(areaName{i},'a','A');
        end
        switch labels_2{i}
            case 'Brainstem'
                areaFile{i} = dir([parentDir.folder, '/BrainstemNavigator/0.9/2a.BrainstemNucleiAtlas_MNI/labels_probabilistic/', ...
                    areaName{i}, '.nii.gz']);
            case 'Diencephalic'
                areaFile{i} = dir([parentDir.folder, '/BrainstemNavigator/0.9/2b.DiencephalicNucleiAtlas_MNI/labels_probabilistic/', ...
                    areaName{i}, '.nii.gz']);
        end
    end
end
