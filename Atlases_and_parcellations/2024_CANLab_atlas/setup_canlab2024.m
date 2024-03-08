% run this script to update *.latest files after modifying the atlas and 
% commit/push your git changes to propogate changes across git clones for 
% all CANLab2023 versions.
%
% This takes on the order of an hour on an i7-12700H (12th gen intel core),
% serially. Parallelization would speed things up but memory overhead is
% substantial, so it's run serially to be more accessible for more users.
% You probably still want a 32G RAM system though.

close all; clear all;

%LIB = '/home/bogdan/.matlab';
LIB = '/dartfs-hpc/rc/home/m/f0042vm/software/';
addpath(genpath([LIB, '/canlab/CanlabCore']))
addpath(genpath([LIB, '/canlab/Neuroimaging_Pattern_Masks']))
%addpath([LIB, '/spm/spm12']);
addpath([LIB, '/spm12']);

for space = {'MNI152NLin6Asym', 'MNI152NLin2009cAsym'}
    for res = [1,2]
        for scale = {'fine','coarse'}
            create_CANLab2024_atlas(space{1},scale{1},res);
        end
    end
end

canlab2023 = load_atlas('canlab2024_coarse_fsl6_2mm');
create_CANLab2023_CIFTI_subctx('MNI152NLin6Asym','coarse',2, canlab2023)
system(which('create_CANLab2024_atlas_cifti.sh'))