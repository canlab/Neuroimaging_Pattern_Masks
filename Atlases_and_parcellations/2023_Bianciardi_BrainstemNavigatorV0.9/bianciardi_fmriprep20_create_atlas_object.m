% if rerun this would benefit from replacing regexp with 
% format_text_letters_only(...,'numbers','cleanup') I think. That's what
% other scripts use, but I didn't know about it when running this the first
% time through.

clear all; close all;

addpath(genpath('~/software/spm12'));

addpath(genpath('~/software/canlab/CanlabCore'))
addpath(genpath('~/software/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('~/software/canlab/MasksPrivate'))

atlas_name = 'bianciardi_fmriprep20';
space_description = 'MNI152NLin2009cAsym';
references = char({'García-Gomar MG, Videnovic A, Singh K, Stauder M, Lewis LD, Wald LL, Rosen BR, Bianciardi M. Disruption of brainstem structural connectivity in RBD using 7 Tesla MRI. Mov Disord. 2021 Dec 29. doi: 10.1002/mds.28895. Online ahead of print. PMID: 34964520',...
    'Singh K, García-Gomar MG, Bianciardi M. Probabilistic Atlas of the Mesencephalic Reticular Formation, Isthmic Reticular Formation, Microcellular Tegmental Nucleus, Ventral Tegmental Area Nucleus Complex, and Caudal-Rostral Linear Raphe Nucleus Complex in Living Humans from 7 Tesla Magnetic Resonance Imaging. Brain Connect. 2021 Oct;11(8):613-623. doi: 10.1089/brain.2020.0975. Epub 2021 biancian 17. PMID: 33926237.',...
    'Singh K, Indovina I, Augustinack JC, Nestor K, García-Gomar MG, Staab JP, Bianciardi M. Probabilistic Template of the Lateral Parabrachial Nucleus, Medial Parabrachial Nucleus, Vestibular Nuclei Complex, and Medullary Viscero-Sensory-Motor Nuclei Complex in Living Humans From 7 Tesla MRI. Front Neurosci. 2020 Jan 23;13:1425. doi: 10.3389/fnins.2019.01425. PMID: 32038134; PMCID: PMC6989551.',...
    'García-Gomar MG, Strong C, Toschi N, Singh K, Rosen BR, Wald LL, Bianciardi M. In vivo Probabilistic Structural Atlas of the Inferior and Superior Colliculi, Medial and Lateral Geniculate Nuclei and Superior Olivary Complex in Humans Based on 7 Tesla MRI. Front Neurosci. 2019 Aug 7;13:764. doi: 10.3389/fnins.2019.00764. PMID: 31440122; PMCID: PMC6694208.',...
    'Bianciardi M, Strong C, Toschi N, Edlow BL, Fischl B, Brown EN, Rosen BR, Wald LL. A probabilistic template of human mesopontine tegmental nuclei from in vivo 7T MRI. Neuroimage. 2018 Apr 15;170:222-230. doi: 10.1016/j.neuroimage.2017.04.070. Epub 2017 May 3. PMID: 28476663; PMCID: PMC5670016.',...
    'Bianciardi M, Toschi N, Edlow BL, Eichner C, Setsompop K, Polimeni JR, Brown EN, Kinney HC, Rosen BR, Wald LL. Toward an In Vivo Neuroimaging Template of Human Brainstem Nuclei of the Ascending Arousal, Autonomic, and Motor Systems. Brain Connect. 2015 Dec;5(10):597-607. doi: 10.1089/brain.2015.0347. Epub 2015 Aug 11. PMID: 26066023; PMCID: PMC4684653.'});

% imort atlas file in MNI152NLin2009cAsym space
% we bianciast use this as a stand in template that we'll modify later, since
% this is the wrong space
bianciaTbl = readtable(which('bianciardi_labels.csv'));
labels = {};
labels_2 = {};
label_descriptions = {};
for i = 1:height(bianciaTbl)
    labels{end+1} = bianciaTbl.left{i};
    label_descriptions{end+1} = bianciaTbl.full_label{i};
    labels_2{end+1} = bianciaTbl.Structure{i};
    if ~isempty(bianciaTbl.right{i})
        label_descriptions{end} = ['Left ',label_descriptions{end}];
        labels{end+1} = bianciaTbl.right{i};
        label_descriptions{end+1} = ['Right ' bianciaTbl.full_label{i}];    
        labels_2{end+1} = bianciaTbl.Structure{i};
    end
end
% fix capitalization
for i = 1:length(label_descriptions), label_descriptions{i}(2:end) = lower(label_descriptions{i}(2:end)); end

% combine data with labels
ref_file = which('fmriprep20_template.nii');
ref = fmri_data(ref_file);
bianciaData = ref;
bianciaAtlas = atlas(bianciaData, ...
    'atlas_name', atlas_name,...
    'labels',labels, ...
    'label_descriptions', label_descriptions, ... 
    'labels_2', labels_2, ...
    'space_description', space_description, ...
    'references',references, 'noverbose');
bianciaAtlas = bianciaAtlas.replace_empty();

% import probability maps
pmap = zeros(size(bianciaData.dat,1), length(bianciaAtlas.labels));
parfor i = 1:length(bianciaAtlas.labels)
    areaName = bianciaAtlas.labels{i};
    switch areaName
        case {'isRt_l', 'isRt_r'}
            areaName = strrep(areaName,'t','T');
        case {'mRt_l'}
            areaName = strrep(areaName,'t','T');
    end
    switch bianciaAtlas.labels_2{i}
        case 'Brainstem'
            file = dir(['BrainstemNavigator/0.9/2a.BrainstemNucleiAtlas_MNI/labels_probabilistic/', ...
                areaName, '_MNI152NLin2009cAsym.nii.gz']);
        case 'Diencephalic'
            file = dir(['BrainstemNavigator/0.9/2b.DiencephalicNucleiAtlas_MNI/labels_probabilistic/', ...
                areaName, '_MNI152NLin2009cAsym.nii.gz']);
    end
    try
        pdata = fmri_data([file.folder, '/', file.name]);   
        system(['rm -f ', file.folder, '/', strrep(file.name,'.gz','')]);
    catch
        error('failed to open file #%d, %s', i, [file.folder, '/', file.name]);
    end

    % super low value errors are likely resampling issues. Not a problem
    % for FSL6 space, but a problem for any other spaces we might resample
    % to in other versions of this script.
    pdata.dat(pdata.dat < 1e-4) = 0;

    pmap(:,i) = pdata.dat;
end
delete(gcp('nocreate'))


% prefix laterality for consistency with other atlases
for i = 1:length(labels), labels{i} = regexprep(labels{i},'(.*)_l$','L_$1'); end
for i = 1:length(labels), labels{i} = regexprep(labels{i},'(.*)_r$','R_$1'); end


bianciaAtlas.probability_maps = sparse(pmap);


% Threshold at probability 0.2 or greater and k = 3 voxels or greater
bianciaAtlas.dat(:) = 0;
bianciaAtlas = threshold(bianciaAtlas, .2, 'k', 3);

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(bianciaAtlas, 'unique', 'overlay', which('fmriprep20_template.nii'));
figure;
% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(bianciaAtlas);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'full2', 'overlay', which('fmriprep20_template.nii'));
    brighten(.6)
    
    o2 = montage(r, o2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

    
end

%% save object

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save([pwd, '/' savename], 'bianciaAtlas');
    
end

%% write - this writes only the label image

if dosave
    
    savename = sprintf('%s_atlas_regions.img', atlas_name);
    bianciaAtlas.fullpath = fullfile(pwd, savename);
    write(bianciaAtlas,'overwrite');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

for i = 1:length(r)
    
    eval([bianciaAtlas.labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save([pwd, '/' savename], 'r', 'region_names', bianciaAtlas.labels{:});

%%
if dosave
    
    figure; han = isosurface(bianciaAtlas);
    
    cellfun(@(x1)set(x1,'FaceAlpha', .5), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1:2), 570, 674]);

    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end
