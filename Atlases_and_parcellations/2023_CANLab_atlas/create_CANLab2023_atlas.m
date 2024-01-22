function atlas_obj = create_CANLab2023_atlas(SPACE, SCALE, res)
% adds bianciardi atlas to canlab2023 scaffolding in the desired space
% Designed for invocation by load_atlas() if atlas with restricted
% bianciardi data isn't already available. This will create it. It depends
% on the bianciardi atlas creation scripts too which will be enchained with
% this if a bianciardi atlas does not yet exist either.
%
% SPACE - MNI152NLin2009cAsym | MNI152NLin6Asym
% SCALE - fine | coarse, describes scale of the parcellation
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

    fprintf('Creating volumetric CANLab2023 %s %s %0.1fmm atlas\n', SCALE, SPACE, res);

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
    atlas_obj = load(sprintf('%s/src/CANLab2023_%s_%s_scaffold.mat',fname.folder,SCALE,SPACE),'canlab');
    atlas_obj = atlas_obj.canlab;

    if res == 2
        template = which(sprintf('%s_T1_2mm.nii.gz',SPACE));
        atlas_obj = atlas_obj.resample_space(template);
        atlas_obj.atlas_name = [atlas_obj.atlas_name '_2mm'];

        biancia = load_atlas(sprintf('bianciardi_%s_%s_2mm',SCALE,alias));
    else
        biancia = load_atlas(sprintf('bianciardi_%s_%s',SCALE,alias));
    end

    % this to handle a bug in the bianciardi files. I'm fixing it in a
    % different matlab session and this line can likely be removed in the
    % future. It's meant to keep this script functional until the fix
    % finishes running
    biancia.label_descriptions = biancia.label_descriptions(:);
    
    biancia.labels_4 = repmat({'Restricted (see BrainstemNavigator v0.9 license)'}, 1, num_regions(biancia));
    biancia.labels_5 = repmat({'Bianciardi brainstem navigator v.0.9'}, 1, num_regions(biancia));
    
    % note that the locus coerulues has more rigorous segmentations based on 
    % T1-TSE sequences, but they produce regions that overlap very well with 
    % the bianciardi atlas' segmentation at this resolution.
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
    error('Deal with the subthalamic atlas. Take bianciardi and merge 1 and 2')
    switch SCALE
        case 'fine'
            biancia_regions = {'PAG', 'SC' 'IC', 'DR', 'MnR', 'SN', 'LC', 'RMg', 'VSM', ...
                'CnF', 'RPa', 'ROb', 'CLi_RLi', 'Rt', 'ION', 'SOC', 'LDTg_CGPn', ...
                'MiTg_PBG', 'PnO_PnC', 'PTg', 'RN', 'SubC', 'Ve', 'STh', 'LPB', 'MPB'};
            biancia = biancia.select_atlas_subset(biancia_regions);
        case 'coarse'
            % we're going to combine related and neighboring areas
            biancia_regions = {'PAG', 'SC' 'IC', 'DR', 'SN', 'VSM', ...
                {'RPa', 'ROb'}, ... % Neighboring medullary raphe nuclei 
                'CLi_RLi', ...
                {'iMRt','PCRtA','sMRt'}, ... % medullary reticular fomrations
                {'isRt','mRt'}, ... mesencephalic regicular formations
                {'ION', 'SOC'}, ... % Inferior and superior olives
                'LDTg_CGPn', 'PnO_PnC', 'MiTg_PBG', 'PTg', 'RN', ...
                {'LC', 'SubC'}, ... % Locus coeruleus and sub coeruleus
                'Ve', 'STh', {'LPB', 'MPB'}, ... % lateral and medial parabrachial nuceli
                'Mn'}; % median and paramedian raphe nuclei
            biancia_names = {'PAG','SC','IC','DR','SN','VSM','RPaOb_B1_B2',...
                'CLi_RLi','MedRt','MesRt','OC','LDTg_CGPn','PnO_PnC',...
                'MiTg_PBG', 'PTg', 'RN', 'LCSubC', 'Ve', 'STh', 'PBN', 'Mn'};
            lateralize_roi = [false,true,true,false,true,true,false,...
                false,true,true,true,true,true,true,true,true,true,true,true,true,...
                false];
            for i = 1:length(biancia_regions)
                ind = find(contains(biancia.labels,biancia_regions{i}));

                if iscell(biancia_regions{i})
                    this_atlas_obj = biancia.select_atlas_subset(biancia_regions{i},'flatten');
                else
                    this_atlas_obj = biancia.select_atlas_subset(biancia_regions(i),'flatten');
                end
                this_atlas_obj.labels = biancia_names(i);

                new_lbls = unique(cellfun(@(x1)regexprep(x1,{'Right ','Left '},''),biancia.label_descriptions(ind),'UniformOutput',false))';
                new_desc = new_lbls{1};
                for j = 2:length(new_lbls)
                    new_desc = [new_desc ', ' new_lbls{j}];
                end
                this_atlas_obj.label_descriptions = {new_desc};
                for lbl = {'labels_2','labels_3','labels_4','labels_5'}
                    if length(biancia.(lbl{1})) == num_regions(biancia)
                        this_atlas_obj.(lbl{1}) = unique(biancia.(lbl{1})(ind));
                    end
                end
                if lateralize_roi(i)
                    this_atlas_obj = lateralize(this_atlas_obj);
                end

                if i == 1
                    biancia_new = this_atlas_obj;
                else
                    biancia_new = biancia_new.merge_atlases(this_atlas_obj);
                end
            end
            biancia = biancia_new;
        otherwise
            error('Did not understand scale %s', SCALE);
    end

    biancia.labels = cellfun(@(x1)(sprintf('Bstem_%s', x1)), biancia.labels, 'UniformOutput', false);
    
    biancia.labels = cellfun(@(x1)(regexprep(x1,'_([LR])_(.*)','_$2_$1')), biancia.labels, 'UniformOutput', false);
    biancia.labels = cellfun(@(x1)(regexprep(x1,'_rh$','_R')), biancia.labels, 'UniformOutput', false);
    biancia.labels = cellfun(@(x1)(regexprep(x1,'_lh$','_L')), biancia.labels, 'UniformOutput', false);
    biancia.labels = cellfun(@(x1)(regexprep(x1,'^([LR])_(.*)','$2_$1')), biancia.labels, 'UniformOutput', false);

    if strcmp(SCALE,'coarse')
        % this atlas needs slack, so we care even about low probability
        % regions, but the brainstem scaffold has p = 0.35, so we increase
        % the pmap here to always exceed the scaffold. By translating and
        % rescaling we preserve the internal rank order of probabilities
        % among bianciardi regions though.
        pmap = biancia.probability_maps;
        rescale = @(x1)(0.351 + x1*(1-0.351));
        biancia.probability_maps(pmap > 0 & pmap <= 0.35) = rescale(pmap(pmap > 0 & pmap <= 0.35));
    end

    atlas_obj = atlas_obj.merge_atlases(biancia);
    atlas_obj.fullpath = '';

    atlas_obj.probability_maps = sparse(double(atlas_obj.probability_maps));
    atlas_obj.references = char(unique(atlas_obj.references));

    timestamp = posixtime(datetime('Now'));
    atlas_obj.additional_info = struct('creation_date', {posixtime(datetime('Now'))});

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
        fprintf('Creating QSIPrep compatable CANLab2023 %s %s %0.1fmm atlas...\n', SCALE, SPACE, res);

        % throw a warning. We don't want to be resamlling atlases again. We
        % should only resample once per atlas and this has already been
        % resampled when its source atlases were projected into whatever
        % space we're operating in now. In particular the Bianciardi atlas
        % is likely to have some misalignment as a result of this second
        % interpolation. We need to integrate sampling to the qsiprep
        % template into the original constituent atlas creation scripts.
        warning('QSIPrep atlas being created needs more work. Resampling non-optimal atlas for now.');
        ref = fmri_data(which('MNI152NLin2009cAsym_1mm_t1s_lps.nii.gz'));
        atlas_obj = atlas_obj.resample_space(ref);
        atlas_obj.fullpath = sprintf('%s/qsiprep/CANLab2023_%s_%s_%dmm_qsiprep.nii', this_dir.folder, SCALE, SPACE, round(res));
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