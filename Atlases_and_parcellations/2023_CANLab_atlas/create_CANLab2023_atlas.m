function atlas_obj = create_CANLab2023_atlas(SPACE, SCALE, res)
% adds bianciardi atlas to canlab2023 scaffolding in the desired space
% Designed for invocation by load_atlas() if atlas with restricted
% bianciardi data isn't already available. This will create it. It depends
% on the bianciardi atlas creation scripts too which will be enchained with
% this if a bianciardi atlas does not yet exist either.
%
% SPACE - MNI152NLin2009cAsym | MNI152NLin6Asym
% SCALE - fine | coarse, describes scale of the parcellation. Generally
%         only affects subcortical structures
% res - 1 | 2. Resolution of final files in mm. If 2 then data is
% additionally thresholded in a number of ways.
    
    switch SPACE
        case 'MNI152NLin2009cAsym'
            alias = 'fmriprep20';
        case 'MNI152NLin6Asym'
            alias = 'fsl6';
        otherwise
            error('No atlas available in %s space', SPACE)
    end

    fprintf('Creating %s volumetric CANLab2023 %s %0.1fmm atlas\n', SCALE, SPACE, res);

    %{
    pmap = fmri_data(which(sprintf('CANLab2023_%s_%s_scaffold.nii.gz',SCALE, SPACE)));
    ref = readtable(which(sprintf('CANLab2023_%s_%s_ref.txt',SCALE, SPACE)));
    lbl = readcell(which(sprintf('CANLab2023_%s_%s_labels.csv',SCALE, SPACE)),'Delimiter',',', 'NumHeaderLines', 1);

    for i = 1:size(lbl,1)
        for j = 1:size(lbl,2)
            if ismissing(lbl{i,j})
                lbl(i,j) = {''};
            end
        end
    end

    atlas_obj = atlas(pmap, ...
        'atlas_name', sprintf('CANLab2023_%s_%s', SCALE, alias), ...
        'labels', lbl(:,1)',...
        'labels_2', lbl(:,3)', ...
        'labels_3', lbl(:,4)', ...
        'labels_4', lbl(:,5)', ...
        'labels_5', lbl(:,6)', ...
        'label_descriptions', lbl(:,2), ...
        'space_description', SPACE,...
        'references', char(ref.references));
    %}
    fname = dir(which('create_CANLab2023_atlas.m'));
    atlas_obj = load(sprintf('%s/src/CANLab2023_%s_scaffold.mat',fname.folder,SPACE),'canlab');
    atlas_obj = atlas_obj.canlab;

    if res == 2
        template = which(sprintf('%s_T1_2mm.nii.gz',SPACE));
        atlas_obj = atlas_obj.resample_space(template);
        atlas_obj.atlas_name = sprintf('%s_%s_%dmm',atlas_obj.atlas_name,SCALE,res);

        biancia = load_atlas(sprintf('bianciardi_%s_2mm',alias));
    else
        atlas_obj.atlas_name = sprintf('%s_%s',atlas_obj.atlas_name,SCALE);
        biancia = load_atlas(sprintf('bianciardi_%s',alias));
    end

    % structs we take from the CIT168 atlas or Kragel2019 (PAG) or LGN/MGN
    % (Morel) that overlap with bianciardi. We prefer Kragel2019 and CIT168
    % to bianciardi because they have open licenses, and we prefer Morel
    % because it has a finer parcellation of the LGN.
    exclude_structs = {'PAG','RN','SN','VTA_PBP','STh', 'LG', 'MG'};
    biancia = biancia.select_atlas_subset(find(~contains(biancia.labels, exclude_structs)));
    
    groupings = {{'DR_B7','MnR_B6_B8','PMnR_B6_B8','CLi_RLi'},...
        {'ROb_B2','RPa_B1','RMg_B3'},...
        {'CnF','IC','SC','MiTg_PBG'},...
        {'isRt','mRta','mRtd','mRtl'},...
        {'PnO_PnC_B5'},...
        {'iMRtl','iMRtm','PCRtA','sMRtl','sMRtm'},...
        {'ION','SOC'},...
        {'LC','SubC'},...
        {'LPB','MPB'},...
        {'LDTg_CGPn', 'PTg'},...
        {'Ve', 'VSM'}};

    labels_3 = {};
    labels_4 = {};
    for i = 1:length(biancia.labels)
        group_ind = cellfun(@(x1)any(contains(biancia.labels{i},x1)),groupings);
        all_group_lbls = biancia.labels(contains(biancia.labels,groupings{group_ind}));
        if contains(biancia.labels{i}, 'L_') && all(cellfun(@(x1)contains(x1,{'L_','R_'}),all_group_lbls))
            side = 'L_';
        elseif contains(biancia.labels{i},'R_') &&  all(cellfun(@(x1)contains(x1,{'L_','R_'}),all_group_lbls))
            side = 'R_';
        else
            side = '';
        end

        % drop L_ or R_ prefixes to allow switch cases to apply bilaterally
        this_label = regexprep(biancia.labels{i},'^[LR]_','');
        switch this_label
            case groupings{1}
                labels_3{end+1} = 'Rostral Raphe (Serotonergic)';
                labels_4{end+1} = 'Midbrain';
            case groupings{2}
                labels_3{end+1} = 'Medullary Raphe (Serotonergic)';
                labels_4{end+1} = 'Medulla';
            case groupings{3}
                labels_3{end+1} = 'Tectum';
                labels_4{end+1} = 'Midbrain';
            case groupings{4}
                labels_3{end+1} = 'Rostral reticular formation';
                labels_4{end+1} = 'Midrain';
            case groupings{5}
                labels_3{end+1} = biancia.labels_2{i};
                labels_4{end+1} = 'Pons';
            case groupings{6}
                labels_3{end+1} = 'Medullary reticular formation';
                labels_4{end+1} = 'Medulla';
            case groupings{7}
                labels_3{end+1} = 'Olivary complex';
                labels_4{end+1} = 'Medulla';
            case groupings{8}
                labels_3{end+1} = 'LC+';
                labels_4{end+1} = 'Pons';
            case groupings{9}
                labels_3{end+1} = 'Parabrachial nuclei';
                labels_4{end+1} = 'Pons';
            case groupings{10}
                labels_3{end+1} = 'Cholinergic nuclei';
                % calling this the pons is a stretch, the PTg is really
                % midbrain, but this way we have the cholinergic nuclei in
                % a single group.
                labels_4{end+1} = 'Pons';
            case groupings{11}
                labels_3{end+1} = 'Cranial nucleu';
                labels_4{end+1} = 'Medullar';
            otherwise
                error('Unexpected Bianciardi area %s', biancia.labels{i});
        end
        % add L_ and R_ prefixes back in
        labels_3{end} = [side, labels_3{end}];
        labels_4{end} = [side, labels_4{end}];
    end
    biancia.labels_3 = labels_3;
    biancia.labels_4 = labels_4;
    biancia.labels_5 = repmat({'Bianciardi brainstem navigator v.0.9'}, 1, num_regions(biancia));
    
    % note that the locus coerulues has more rigorous segmentations based on 
    % T1-TSE sequences, but they produce regions that overlap very well with 
    % the bianciardi atlas' segmentation at this resolution.
    % Re: updates from CANLab2018,
    % MnR replaces NCS
    % VSM replaces dmnx_nts
    % RMg replaces ncs_B6_B8
    % NRP (nucleus raphe pontis) - dropped in CANLab2023
    % nuc_ambiguus - dropped in CANLab2023
    % medullary_raphe - replaced by RPa and ROb
    % spinal_trigeminal - temporarily excluded. Can be added later with
    %   transformation from SUIT brainstem space (transforms into MNI available
    %   from Diedrichsen lab)
    % CLi_RLi - new in CANLab2023
    % iMRt, isRt, mRt, PCRtA, PnO_PnC, sMRt - reticular formations, new in CANLab2023
    % SOC, ION - superior and inferior olives, new in CANLab2023
    % MiTg_PBG
    % PDT_CGPn
    % PTg
    % RN - red nucleus new in CANLab2023
    % SubC - subcoeruleus new in CANLab2023
    % Ve - vestibular nucle icomplex
    % STh - Subthalamic nucelus new in CANLab2023
    biancia.labels = cellfun(@(x1)(sprintf('Bstem_%s', x1)), biancia.labels, 'UniformOutput', false);
    
    biancia.labels = cellfun(@(x1)(regexprep(x1,'_([LR])_(.*)','_$2_$1')), biancia.labels, 'UniformOutput', false);
    biancia.labels = cellfun(@(x1)(regexprep(x1,'_rh$','_R')), biancia.labels, 'UniformOutput', false);
    biancia.labels = cellfun(@(x1)(regexprep(x1,'_lh$','_L')), biancia.labels, 'UniformOutput', false);
    biancia.labels = cellfun(@(x1)(regexprep(x1,'^([LR])_(.*)','$2_$1')), biancia.labels, 'UniformOutput', false);

    % merge biancia with other relevant structures
    biancia = atlas_obj.select_atlas_subset(exclude_structs).merge_atlases(biancia);
    % save complementary atlas_obj
    atlas_obj = atlas_obj.select_atlas_subset(find(~contains(atlas_obj.labels, exclude_structs)));
    % renormalize
    total_p = sum(biancia.probability_maps,2);
    renorm = total_p > 1;
    biancia.probability_maps(renorm,:) = biancia.probability_maps(renorm,:)./total_p(renorm);

    % merge augmented biancia with the rest of the atlas, replacing
    % anything that overlaps (shen regions primarily)
    atlas_obj = biancia.merge_atlases(atlas_obj,'noreplace');
    atlas_obj.fullpath = '';

    atlas_obj.probability_maps = sparse(double(atlas_obj.probability_maps));
    atlas_obj.references = char(unique(atlas_obj.references,'rows'));

    timestamp = posixtime(datetime('Now'));
    atlas_obj.additional_info = struct('creation_date', {posixtime(datetime('Now'))});

    if strcmp(SCALE,'coarse')
        atlas_obj = atlas_obj.downsample_parcellation('labels_2');
    end

    this_dir = dir(which('create_CANLab2023_atlas.m'));
    savename = sprintf('%s_atlas_object.mat', atlas_obj.atlas_name);
    save([this_dir.folder, '/' savename], 'atlas_obj');

    % we can't upload the mat file to github due to licensing issues, but
    % we can upload a timestamp that will flag out of date versions and
    % cause other uesrs to recreate the atlas when appropriate.
    fid = fopen(sprintf('%s/%s_atlas_object.latest', this_dir.folder, atlas_obj.atlas_name),'w+');
    fprintf(fid,'%f',timestamp);
    fclose(fid);

    if any(ismember(SPACE,{'MNI152NLin6Asym'})) && strcmp(SCALE,'coarse') && res == 2
        % save subcortical volume for use in CIFTI creation
        fprintf('Creating CIFTI atlas...\n')
        create_CANLab2023_CIFTI_subctx(SPACE,SCALE,res,atlas_obj);
    end

    if any(ismember(SPACE,{'MNI152NLin2009cAsym'})) && strcmp(SCALE,'coarse') && res == 1
        fprintf('Creating QSIPrep compatable CANLab2023 %s %s %0.1fmm atlas...\n', SPACE, SCALE, res);

        % throw a warning. We don't want to be resampling atlases again. We
        % should only resample once per atlas and this has already been
        % resampled when its source atlases were projected into whatever
        % space we're operating in now. In particular the Bianciardi atlas
        % is likely to have some misalignment as a result of this second
        % interpolation. We need to integrate sampling to the qsiprep
        % template into the original constituent atlas creation scripts.
        warning('QSIPrep atlas being created needs more work. Resampling non-optimal atlas for now.');
        ref = fmri_data(which('MNI152NLin2009cAsym_1mm_t1s_lps.nii.gz'));
        atlas_obj = atlas_obj.resample_space(ref).threshold(0.2);
        atlas_obj.fullpath = sprintf('%s/qsiprep/CANLab2023_%s_%s_%dmm_qsiprep.nii', this_dir.folder, SPACE, SCALE, round(res));
        atlas_obj.write('overwrite');
        
        % zero sform for compliance with qsiprep reqs. See here under 
        % "Using custom atlases",
        % https://qsiprep.readthedocs.io/en/latest/reconstruction.html
        V = niftiread(atlas_obj.fullpath);
        info = niftiinfo(atlas_obj.fullpath);
        info.TransformName = 'Qform';
        niftiwrite(V,atlas_obj.fullpath,info);
        
        gzip(atlas_obj.fullpath);
        delete(atlas_obj.fullpath);

        % create atlas_config.json file
        jsonstruct = struct(sprintf('canlab2023_%s',SCALE),...
            struct('file', {[atlas_obj.fullpath, '.gz']}, ...
                'node_names', {atlas_obj.labels}, ...
                'node_ids', {1:length(atlas_obj.labels)}));
        json_file = sprintf('%s/qsiprep/atlas_config.json', this_dir.folder);
        if exist(json_file,'file')
            fid = fopen(json_file);
            oldtxt = fread(fid,inf);
            fclose(fid);

            jsonOld = jsondecode(char(oldtxt'));
            old_fields = fieldnames(jsonOld);
            new_fields = fieldnames(jsonstruct);
            if any(ismember(old_fields,fieldnames(jsonstruct)))
                warning('Will overwrite old JSON atlas entry %s\n', old_fields{ismember(old_fields, new_fields)});
                warning('Copying old json from %s to %s',json_file, [json_file '_bak']);

                copyfile(json_file, [json_file '_bak']);
            end
            for fname = old_fields
                jsonstruct.(fname{1}) = jsonOld.(fname{1});
            end
        end
        jsontxt = jsonencode(jsonstruct, 'PrettyPrint', true);

        fid = fopen(json_file,'w+');
        fprintf(fid, '%s', jsontxt);
        fclose(fid);

        fprintf('Wrote qsiprep compatable atlas to %s/qsiprep/\n', this_dir.folder);
    end
end