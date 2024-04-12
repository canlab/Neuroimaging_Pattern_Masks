% run this script to update *.latest files after modifying the atlas and 
% commit/push your git changes to propogate changes across git clones for 
% all CANLab2023 versions.
%
% This takes on the order of an hour on an i7-12700H (12th gen intel core),
% serially. Parallelization would speed things up but memory overhead is
% substantial, so it's run serially to be more accessible for more users.
% You probably still want a 32G RAM system though.

close all; clear all;

LIB = '/home/bogdan/.matlab';
%LIB = '/dartfs-hpc/rc/home/m/f0042vm/software/';
addpath(genpath([LIB, '/canlab/CanlabCore']))
addpath(genpath([LIB, '/canlab/Neuroimaging_Pattern_Masks']))
addpath([LIB, '/spm/spm12']);
%addpath([LIB, '/spm12']);

system(which('create_CANLab2024_atlas_prep.sh'));
for space = {'MNI152NLin6Asym', 'MNI152NLin2009cAsym'}
    for res = [1,2]
        for scale = {'fine','coarse'}
            create_CANLab2024_atlas(space{1},scale{1},res);
        end
    end
end

canlab2024 = load_atlas('canlab2024_coarse_fsl6_2mm');
create_CANLab2024_CIFTI_subctx('MNI152NLin6Asym','coarse',2, canlab2024)
system(which('create_CANLab2024_atlas_cifti.sh'))

%% save label files
canlab2024 = load_atlas('canlab2024_fine');
lbl_tbl = table(canlab2024.labels', canlab2024.labels_2', canlab2024.labels_3', ...
    canlab2024.labels_4', canlab2024.labels_5', canlab2024.label_descriptions,...
    'VariableNames',{'labels_1','labels_2','labels_3','labels_4','source','label_description'});

outdir = dir(which('evaluate_canlab2024.m'));
writetable(lbl_tbl,fullfile(outdir.folder,'canlab2024_labels.csv'),'delimiter',',')

%% save references
outdir = dir(which('evaluate_canlab2024.m'));
fid = fopen(fullfile(outdir.folder,'canlab2024_references.txt'),'w+');
for i = 1:size(canlab2024.references,1)
    fprintf(fid,'%s\n',canlab2024.references(i,:));
end
fclose(fid);

%% save probability maps
% this requires a ton of memory to run and the result cannot be redistributed
% anyway, only uncomment if you need nifti files for personal use.
%{
outdir = dir(which('evaluate_canlab2024.m'));
for res = [2,1]
    for alias = {'fmriprep20', 'fsl6'}
        this_atlas = load_atlas(sprintf('canlab2024_fine_%s_%dmm',alias{1},res));
        pmap = fmri_data(this_atlas);
        pmap.dat = single(full(this_atlas.probability_maps));
        switch alias{1}
            case 'fmriprep20'
                space = 'MNI152NLin2009cAsym';
            case 'fsl6'
                space = 'MNI152NLin6Asym';
            otherwise
                error('Unrecognized alias');
        end
        pmap.fullpath = fullfile(outdir.folder,sprintf('CANLab2024_%s_%dmm.nii',space,res));
        pmap.write('overwrite');
        gzip(pmap.fullpath);
        delete(pmap.fullpath);
    end
end
%}