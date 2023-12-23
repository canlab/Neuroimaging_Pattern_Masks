close all; clear all;

addpath('/dartfs-hpc/rc/home/m/f0042vm/software/spm12')
addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/canlab/CanlabCore/'))

SPACE = 'MNI152NLin2009cAsym';

lh_lbls = readtable('lctx_labels.txt');
rh_lbls = readtable('rctx_labels.txt');
lbls = [lh_lbls{:,1}; rh_lbls{:,1}];

pmap = {};
for study = {'bmrk5', 'paingen', 'spacetop'}
    for hemi = {'lh' 'rh'}
        pdata = fmri_data(sprintf('%s_%s_%s.nii.gz', hemi{1}, study{1}, SPACE));
        
        pmap{end+1} = zeros(size(pdata.dat,1), length(lbls)/2);
        for i = 1:length(lbls)/2
            pmap{end}(:,i) = mean(pdata.dat == i,2);
        end 
    end
    % combine left and right hemispheres
    pmap{end-1} = cat(2,pmap{end-1:end});
    pmap(end) = [];
end

pdata = pdata.get_wh_image(1);

glasser = cell(length(pmap),1);
for i = 1:length(pmap)
    pdata.dat = pmap{i};
    glasser{i} = atlas(pdata);
    glasser{i}.labels = lbls';
    glasser{i}.probability_maps = sparse(glasser{i}.probability_maps);
end

% save individual atlases for inspection
glasser{1}.fullpath = sprintf('bmrk5_%s_atlas.nii', SPACE);
glasser{2}.fullpath = sprintf('paingen_%s_atlas.nii', SPACE);
glasser{3}.fullpath = sprintf('spacetop_%s_atlas.nii', SPACE);
glasser{1}.threshold(0.2).write('overwrite');
glasser{2}.threshold(0.2).write('overwrite');
glasser{3}.threshold(0.2).write('overwrite');
gzip(glasser{1}.fullpath)
gzip(glasser{2}.fullpath)
gzip(glasser{3}.fullpath)

% save subset atlases for evaluation on held out data
glasser_joint = cell(length(pmap),2);
for i = 1:length(pmap)
    these_ind = find(~ismember(1:length(pmap),i));
    pdata.dat = sparse(mean(cat(3,pmap{these_ind}),3));
    glasser_joint{i} = atlas(pdata);
    glasser_joint{i}.labels = lbls';
end

glasser_joint{1}.fullpath = sprintf('no_bmrk5_%s_atlas.nii', SPACE);
glasser_joint{2}.fullpath = sprintf('no_paingen_%s_atlas.nii', SPACE);
glasser_joint{3}.fullpath = sprintf('no_spacetop_%s_atlas.nii', SPACE);
glasser_joint{1}.threshold(0.2).write('overwrite');
glasser_joint{2}.threshold(0.2).write('overwrite');
glasser_joint{3}.threshold(0.2).write('overwrite');
gzip(glasser_joint{1}.fullpath)
gzip(glasser_joint{2}.fullpath)
gzip(glasser_joint{3}.fullpath)

%% generate final mean atlas
labels = lbls;
labels = format_text_letters_only(labels, 'numbers', 'cleanup'); % Replace chars we don't want
labels = strrep(labels, '_ROI', '');

labels = regexprep(labels, '(^[LR])_(.*)', 'Ctx_$2_$1');

lbl_descrip = table2cell(readtable('src/label_descriptions.csv','Delimiter','\t'));
lbl_descrip = [cellfun(@(x1)([x1, ' (Left)']), lbl_descrip, 'UniformOutput', false);  cellfun(@(x1)([x1, ' (Right)']), lbl_descrip, 'UniformOutput', false)];

references = char([{'Glasser, Matthew F., Timothy S. Coalson, Emma C. Robinson, Carl D. Hacker, John Harwell, Essa Yacoub, Kamil Ugurbil, et al. 2016. A Multi-Modal Parcellation of Human Cerebral Cortex. Nature 536 (7615): 171?78.'; ...
    'Wu J, Ngo GH, Greve D, Li J, He T, Fischl B, Eickhoff SB, Yeo T. Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems. 2018. Human Brain Mapping 39(9) 3793-3808. DOI: 10.1002/hbm.24213'}]);

meanPmap = mean(cat(3,pmap{:}),3);
pdata.dat = meanPmap;
atlas_obj = atlas(pdata, ...
    'atlas_name', sprintf('glasser_%s',SPACE), ...
    'labels', labels',...
    'label_descriptions', lbl_descrip,...
    'space_description', SPACE, ...
    'references', references, 'noverbose');
atlas_obj.probability_maps = sparse(atlas_obj.probability_maps);


savename = sprintf('%s_atlas_object.mat', atlas_obj.atlas_name);
save(savename, 'atlas_obj');


% save nifti version
nii = fmri_data(atlas_obj);
nii.fullpath = sprintf('glasser_%s_atlas.nii', SPACE);
nii.dat = single(full(atlas_obj.probability_maps));
nii.write('overwrite');
gzip(nii.fullpath);
delete(nii.fullpath);
fid = fopen(sprintf('glasser_%s_atlas_labels.txt',SPACE),'w+');
fprintf(fid,'%s\n', atlas_obj.labels{:})
fclose(fid);

% save thresholded convenience labels
nii = fmri_data(atlas_obj.threshold(0.2));
nii.fullpath = sprintf('glasser_%s_atlas_p20.nii', SPACE);
nii.write('overwrite');
gzip(nii.fullpath);