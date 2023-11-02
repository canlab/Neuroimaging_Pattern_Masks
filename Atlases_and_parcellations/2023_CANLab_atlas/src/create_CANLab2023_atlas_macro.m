clear all; close all;

LIB = '/dartfs-hpc/rc/home/m/f0042vm/software';
%LIB = '/home/bogdan/.matlab';
ROOT = [LIB, '/canlab/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2023_CANLab_atlas/'];

addpath([LIB, '/spm12']);
%addpath([LIB, '/spm/spm12']);

addpath(genpath([LIB, '/canlab/CanlabCore']))
addpath(genpath([LIB, '/canlab/Neuroimaging_Pattern_Masks']))
addpath(genpath([LIB, '/canlab/MasksPrivate']))

%for SPACE = {'MNI152NLin6Asym', 'MNI152NLin2009cAsym'}
for SPACE = {'MNI152NLin2009cAsym'}
    SPACE = SPACE{1};
    for SCALE = {'fine','coarse'}
        SCALE = SCALE{1};
        create_CANLab2023_unrestricted
        switch SCALE
            case 'fine'
                fine = true;
            case 'coarse'
                fine = false;
        end
        for res = 1:2
            if res == 1
                bianciardi_create_atlas_obj(SPACE,fine);
            elseif res == 2
                bianciardi_create_atlas_obj([SPACE,'_2mm'],fine);
            end
            create_CANLab2023_atlas(SPACE,SCALE,res)
        end
    end
end