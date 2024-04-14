function atlas_obj = create_CANLab2024_atlas(SPACE, SCALE, res)
% adds bianciardi atlas to canlab2024 scaffolding in the desired space
% Designed for invocation by load_atlas() if atlas with restricted
% bianciardi data isn't already available. This will create it. It depends
% on the bianciardi atlas creation scripts too which will be enchained with
% this if a bianciardi atlas does not yet exist either.
%
% SPACE - MNI152NLin2009cAsym | MNI152NLin6Asym
% SCALE - fine | coarse, describes scale of the parcellation. Generally
%         only affects subcortical structures
% res - 1 | 2. Resolution of final files in mm. If 2 then data is
% versionUpdate - true/false. Set to true if you've just modified the atlas
%   somehow. If this script is invoked simply to build a copy of the atlas
%   that you don't have locally (e.g. from load_atlas()), or to pull an 
% additionally thresholded in a number of ways.
    

% ToDo MnR esems to be lost in the fine 2mm MNI152NLin6Asym atlas. Might
% make sense to subsume this into PMnR for all 2mm atlases

    switch SPACE
        case 'MNI152NLin2009cAsym'
            alias = 'fmriprep20';
        case 'MNI152NLin6Asym'
            alias = 'fsl6';
        otherwise
            error('No atlas available in %s space', SPACE)
    end

    fprintf('Creating %s volumetric CANLab2024 %s %0.1fmm atlas\n', SCALE, SPACE, res);

    %{
    pmap = fmri_data(which(sprintf('CANLab2024_%s_%s_scaffold.nii.gz',SCALE, SPACE)));
    ref = readtable(which(sprintf('CANLab2024_%s_%s_ref.txt',SCALE, SPACE)));
    lbl = readcell(which(sprintf('CANLab2024_%s_%s_labels.csv',SCALE, SPACE)),'Delimiter',',', 'NumHeaderLines', 1);

    for i = 1:size(lbl,1)
        for j = 1:size(lbl,2)
            if ismissing(lbl{i,j})
                lbl(i,j) = {''};
            end
        end
    end

    atlas_obj = atlas(pmap, ...
        'atlas_name', sprintf('CANLab2024_%s_%s', SCALE, alias), ...
        'labels', lbl(:,1)',...
        'labels_2', lbl(:,3)', ...
        'labels_3', lbl(:,4)', ...
        'labels_4', lbl(:,5)', ...
        'labels_5', lbl(:,6)', ...
        'label_descriptions', lbl(:,2), ...
        'space_description', SPACE,...
        'references', char(ref.references));
    %}
    atlas_obj = load_atlas(sprintf('openCANLab2024_%s_%s_%dmm', SCALE, alias, res));
    atlas_obj.atlas_name = strrep(atlas_obj.atlas_name, 'open', '');

    % drop harvard AAN regions
    switch SCALE
        case 'fine'
            source_lbls = 'labels_5';
            % this should correspond to labels_2 of the subcoreulus
            atlas_obj.labels_2(contains(atlas_obj.labels,'LC_L')) = {strrep(atlas_obj.labels_2{contains(atlas_obj.labels,'LC_L')},'LC_L','LC+_L')};
            atlas_obj.labels_2(contains(atlas_obj.labels,'LC_R')) = {strrep(atlas_obj.labels_2{contains(atlas_obj.labels,'LC_R')},'LC_R','LC+_R')};
        case 'coarse'
            source_lbls = 'labels_4';
            % this should correspond to labels of the subcoreulus
            atlas_obj.labels(contains(atlas_obj.labels,'LC_L')) = {strrep(atlas_obj.labels{contains(atlas_obj.labels,'LC_L')},'LC_L','LC+_L')};
            atlas_obj.labels(contains(atlas_obj.labels,'LC_R')) = {strrep(atlas_obj.labels{contains(atlas_obj.labels,'LC_R')},'LC_R','LC+_R')};
    end
    atlas_obj = atlas_obj.select_atlas_subset(find(~contains(atlas_obj.(source_lbls),'Harvard')));

    % drop Levinson-Bari NTS
    atlas_obj = atlas_obj.select_atlas_subset(find(~contains(atlas_obj.labels,'NTS')));

    if res == 2
        biancia = load_atlas(sprintf('bianciardi_%s_2mm',alias));
    else
        biancia = load_atlas(sprintf('bianciardi_%s',alias));
    end

    biancia = biancia.resample_space(atlas_obj);

    % structs we take from the CIT168 atlas, Kragel2019 (PAG), Iglesias2018 
    % (LGN/MGN) or Levinson-Bari (DR, LC and VTS analog NTS)  that overlap 
    % with bianciardi. We prefer these to bianciardi because they have open 
    % licenses.
    % We might want to incorporate the Morel LGN into the Iglesias
    % probabilities down the line if we can figure out some clever way to
    % create hybrid regions (borrowing probabilities from one with
    % subdivisions from the other), but for now we exclude this.
    % note, keep NTS here because we can use this both ways, both with
    % biancia and with atlas_obj
    exclude_structs = {'PAG', 'RN', 'SN', 'VTA_PBP', 'STh', 'LG', 'MG', 'DR_B7', 'LC', 'NTS'};
    biancia = biancia.select_atlas_subset(find(~contains(biancia.labels, exclude_structs)));
    biancia_refs = biancia.labels_4;
    
    % diagnostic code to plot the overlap of biancia with substitute atlases
    %{
    cit = atlas_obj.select_atlas_subset({'RN','SN','VTA','STh'});
    kragel = atlas_obj.select_atlas_subset({'PAG'});
    iglesias = atlas_obj.select_atlas_subset({'LGN','MGN'});
    levinsonbari = atlas_obj.select_atlas_subset({'DR', 'LC', 'NTS'});
    substitute_regions = [cit, kragel, iglesias, levinsonbari];
    
    if res == 2
        bianciaFull = load_atlas(sprintf('bianciardi_%s_2mm',alias));
    else
        bianciaFull = load_atlas(sprintf('bianciardi_%s',alias));
    end
    biancia1 = bianciaFull.select_atlas_subset(find(contains(bianciaFull.labels, exclude_structs)));
    biancia2 = bianciaFull.select_atlas_subset(find(~contains(bianciaFull.labels, exclude_structs)));
    bianciaFull = [biancia1, biancia2];
    %}

    groupings = {{'DR_B7','MnR_B6_B8','PMnR_B6_B8','CLi_RLi'},...
        {'ROb_B2','RPa_B1','RMg_B3'},...
        {'IC','SC','MiTg_PBG'},...
        {'isRt','mRta','mRtd','mRtl'},...
        {'PnO_PnC_B5'},...
        {'iMRtl','iMRtm','PCRtA','sMRtl','sMRtm'},...
        {'ION','SOC'},...
        {'LC','SubC'},...
        {'LPB','MPB'},...
        {'LDTg_CGPn', 'PTg'},...
        {'Ve', 'VSM','NTS'},... % we may or may not include NTS, but it doesn't hurt to account for it here even if it's missing
        {'CnF'}};

    [labels_3, labels_4, labels_5] = deal({});
    for i = 1:length(biancia.labels)
        group_ind = cellfun(@(x1)any(contains(biancia.labels{i},x1)),groupings);
        all_group_lbls = biancia.labels(contains(biancia.labels,groupings{group_ind}));
        if contains(biancia.labels{i}, 'L_') && all(cellfun(@(x1)contains(x1,{'L_','R_'}),all_group_lbls))
            side = '_L';
        elseif contains(biancia.labels{i},'R_') &&  all(cellfun(@(x1)contains(x1,{'L_','R_'}),all_group_lbls))
            side = '_R';
        else
            side = '';
        end

        % drop L_ or R_ prefixes to allow switch cases to apply bilaterally
        this_label = regexprep(biancia.labels{i},'^[LR]_','');
        switch this_label
            case groupings{1}
                % this should match DR labels added to Levinson-Bari
                % regions in create_brainstem2024_atlas_unrestricted.m
                labels_3{end+1} = 'Rostral Raphe (Serotonergic)';
                labels_4{end+1} = 'Midbrain';
            case groupings{2}
                biancia.labels_2{i} = 'RObPaMg';
                labels_3{end+1} = 'Medullary Raphe (Serotonergic)';
                labels_4{end+1} = 'Medulla';
            case groupings{3}
                labels_3{end+1} = 'Tectum';
                labels_4{end+1} = 'Midbrain';
            case groupings{4}
                labels_3{end+1} = 'Rostral reticular formation';
                labels_4{end+1} = 'Midbrain';
            case groupings{5}
                labels_3{end+1} = biancia.labels_2{i};
                labels_4{end+1} = 'Pons';
            case groupings{6}
                labels_3{end+1} = 'Medullary reticular formation';
                labels_4{end+1} = 'Medulla';
            case groupings{7}
                % SOC is too small at level 2, so let's merge it
                biancia.labels_2{i} = strrep(biancia.labels_2{i},'ION','OC');
                biancia.labels_2{i} = strrep(biancia.labels_2{i},'SOC','OC');
                labels_3{end+1} = 'Olivary complex';
                labels_4{end+1} = 'Medulla';
            case groupings{8}
                % this should match LC labels added to Levinson-Bari
                % regions in create_brainstem2024_atlas_unrestricted.m
                biancia.labels_2{i} = strrep(biancia.labels_2{i},'SubC','LC+');
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
                labels_3{end+1} = 'Cranial_nuclei';
                labels_4{end+1} = 'Medulla';
            case groupings{12}
                biancia.labels_2{i} = 'PAG';
                labels_3{end+1} = 'PAG';
                labels_4{end+1} = 'Midbrain';
                side = ''; % delateralize the cueniform due to merger with PAG
            otherwise
                error('Unexpected Bianciardi area %s', this_label);
        end
        % add L_ and R_ prefixes back in
        labels_3{end} = [labels_3{end}, side];
        labels_4{end} = [labels_4{end}, side];
        labels_5{end+1} = ['Bianciardi brainstem navigator v.0.9 (ref: ', biancia_refs{i}, ')'];
    end
    biancia.labels_3 = labels_3;
    biancia.labels_4 = labels_4;
    biancia.labels_5 = labels_5;
    
    % note that the locus coerulues has more rigorous segmentations based on 
    % T1-TSE sequences, but they produce regions that overlap very well with 
    % the bianciardi atlas' segmentation at this resolution.
    % Re: updates from CANLab2018,
    % MnR replaces NCS
    % NTS replaces dmnx_nts (VSM in CANLab2023)
    % RMg replaces ncs_B6_B8
    % NRP (nucleus raphe pontis) - dropped in CANLab2023/24
    % nuc_ambiguus - dropped in CANLab2023/24
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
    biancia.labels = cellfun(@(x1)(sprintf('BStem_%s', x1)), biancia.labels, 'UniformOutput', false);
    biancia.labels_2 = cellfun(@(x1)(sprintf('BStem_%s', x1)), biancia.labels_2, 'UniformOutput', false);
    for lbl = {'labels','labels_2'}        
        biancia.(lbl{1}) = cellfun(@(x1)(regexprep(x1,'_([LR])_(.*)','_$2_$1')), biancia.(lbl{1}), 'UniformOutput', false);
        biancia.(lbl{1}) = cellfun(@(x1)(regexprep(x1,'_rh$','_R')), biancia.(lbl{1}), 'UniformOutput', false);
        biancia.(lbl{1}) = cellfun(@(x1)(regexprep(x1,'_lh$','_L')), biancia.(lbl{1}), 'UniformOutput', false);
        biancia.(lbl{1}) = cellfun(@(x1)(regexprep(x1,'^([LR])_(.*)','$2_$1')), biancia.(lbl{1}), 'UniformOutput', false);
    end

    % merge biancia with other relevant structures
    other_bstem = atlas_obj.select_atlas_subset(exclude_structs);
    % save complementary atlas_obj
    atlas_obj = atlas_obj.select_atlas_subset(find(~contains(atlas_obj.labels, exclude_structs)));
    % renormalize
    total_p = sum(biancia.probability_maps,2) + sum(other_bstem.probability_maps,2);
    renorm = total_p > 1;
    biancia.probability_maps(renorm,:) = biancia.probability_maps(renorm,:)./total_p(renorm);
    other_bstem.probability_maps(renorm,:) = other_bstem.probability_maps(renorm,:)./total_p(renorm);

    if strcmp(SCALE,'coarse')
        biancia = biancia.downsample_parcellation('labels_2');
        
        % merge subcoreulus into coreulus, while keeping indexing unchanged
        LC_ind = find(contains(other_bstem.labels,{'LC+_R'}));
        other_bstem.probability_maps(:,LC_ind) = other_bstem.probability_maps(:,LC_ind) + biancia.select_atlas_subset({'LC+_R'}).probability_maps;
        other_bstem.label_descriptions(LC_ind) = {'Locus coeruleus and subcoeruleus (right)'};

        LC_ind = find(contains(other_bstem.labels,{'LC+_L'}));
        other_bstem.probability_maps(:,LC_ind) = other_bstem.probability_maps(:,LC_ind) + biancia.select_atlas_subset({'LC+_L'}).probability_maps;
        other_bstem.label_descriptions(LC_ind) = {'Locus coeruleus and subcoeruleus (left)'};

        other_bstem = other_bstem.probability_maps_to_region_index;

        biancia = biancia.select_atlas_subset(find(~contains(biancia.labels,'LC+')));
    end

    
    biancia = [other_bstem, biancia];

    %% adjust shen regions so they're always less than biancia regions
    % (since the shen brainstem regions are basically fillers
    shen_regions = contains(atlas_obj.(source_lbls),'Shen');
    % extract shen_regions
    shen = atlas_obj.select_atlas_subset(find(shen_regions));
    atlas_obj = atlas_obj.select_atlas_subset(find(~shen_regions));
    % reconcile biancia regions with nonshen part of atlas
    atlas_obj = atlas_obj.merge_atlases(biancia);
    % renorm
    total_p = sum(atlas_obj.probability_maps,2);
    renorm = total_p > 1;
    atlas_obj.probability_maps(renorm,:) = atlas_obj.probability_maps(renorm,:)./total_p(renorm);

    % renormalize shen to sum to whatever probability density is left over
    % after accounting for other regions
    resid_p = 1 - sum(atlas_obj.probability_maps,2);
    total_p = sum(shen.probability_maps,2);
    new_total_p = min([sum(shen.probability_maps,2), resid_p],[],2);
    renorm = total_p > new_total_p;
    s = new_total_p(renorm)./total_p(renorm);
    s(s<0 & s<-1e-7) = 0; % floating point errors
    shen.probability_maps(renorm,:) = shen.probability_maps(renorm,:).*s;
    
    % add shen back
    atlas_obj = atlas_obj.merge_atlases(shen);

    %% save data
    atlas_obj.fullpath = '';

    atlas_obj.probability_maps = sparse(double(atlas_obj.probability_maps));

    if strcmp('MNI152NLin6Asym',SPACE) && res == 2 && strcmp(SCALE,'fine')
        % hacky fix for LC_L region that does not survive at this
        % resolution in this atlas. At coarse scale we don't have LC
        % anymore, and instead combine LC and SubC into LC+, so we don't
        % encounter this issue.

        canlab_LC_ind = contains(atlas_obj.labels,{'LC_L'});
    
        max_ind = find(max(atlas_obj.probability_maps(:,canlab_LC_ind)) == atlas_obj.probability_maps(:,canlab_LC_ind));
        % make sure we haven't resampled in some weird way and this area 
        % still has probability assigned to LC after merging the bianciardi 
        % atlas into our atlas_obj
        assert(atlas_obj.probability_maps(max_ind,canlab_LC_ind) ~= 0);
        
        LC_pval = full(atlas_obj.probability_maps(max_ind,canlab_LC_ind));
        canlab_not_LC = find(~canlab_LC_ind);
        other_pmaps = full(atlas_obj.probability_maps(max_ind,canlab_not_LC));
        assert(sum(other_pmaps > LC_pval) == 1) % make sure only one region exceeds LC value
        % the above region should be a shen region
        bad_region = canlab_not_LC(other_pmaps > LC_pval);
        
        % find the cumulative probability of the LC and "bad" region
        total_p = sum(atlas_obj.probability_maps(max_ind,[bad_region, find(canlab_LC_ind)]));        
        
        % modify probability to something that won't overrule LC, splitting
        % the difference just marginally.
        new_p = total_p/2;
        atlas_obj.probability_maps(max_ind,bad_region) = new_p-0.01;
        atlas_obj.probability_maps(max_ind,find(canlab_LC_ind)) = new_p+0.01;

        % regenerate pmaps
        atlas_obj = atlas_obj.probability_maps_to_region_index();
    end

    atlas_obj.references = char(unique(atlas_obj.references,'rows'));

    hash = DataHash(atlas_obj);
    atlas_obj.additional_info = struct('creation_date', {posixtime(datetime('Now'))},...
        'hash',{hash});

    [~,~,~,missing_regions] = atlas_obj.check_properties();
    if ~isempty(missing_regions)
        error('Some regions are missing from the atlas file');
    end

    this_dir = dir(which('create_CANLab2024_atlas.m'));
    savename = sprintf('%s_atlas_object.mat', atlas_obj.atlas_name);
    save([this_dir.folder, '/' savename], 'atlas_obj');

    % we can't upload the mat file to github due to licensing issues, but
    % we can upload a timestamp that will flag out of date versions and
    % cause other uesrs to recreate the atlas when appropriate.
    fid = fopen(sprintf('%s/%s_atlas_object.latest', this_dir.folder, atlas_obj.atlas_name),'w+');
    fprintf(fid,'%s',hash);
    fclose(fid);

    if any(ismember(SPACE,{'MNI152NLin6Asym'})) && strcmp(SCALE,'coarse') && res == 2
        % save subcortical volume for use in CIFTI creation
        fprintf('Creating CIFTI atlas...\n')
        create_CANLab2024_CIFTI_subctx(SPACE,SCALE,res,atlas_obj);
    end

    if any(ismember(SPACE,{'MNI152NLin2009cAsym'})) && strcmp(SCALE,'coarse') && res == 1
        fprintf('Creating QSIPrep compatable CANLab2024 %s %s %0.1fmm atlas...\n', SPACE, SCALE, res);

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
        atlas_obj.fullpath = sprintf('%s/qsiprep/CANLab2024_%s_%s_%dmm_qsiprep.nii', this_dir.folder, SPACE, SCALE, round(res));
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
        jsonstruct = struct(sprintf('canlab2024_%s',SCALE),...
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
