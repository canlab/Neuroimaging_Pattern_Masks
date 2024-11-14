close all; clear all;

addpath('/home/bogdan/.matlab/spm/spm12')
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore/'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks/'))

space_description = 'MNI152NLin6Asym';
alias = 'fsl6';
atlas_name = sprintf('destrieux_%s',alias);
references = char({'Fischl B, var der Kouwe A, Destrieux C, Halgren E, Segonne F, Salat DH, Busa E, Seidman LJ, Goldstein J, Kennedy D, Caviness V, Markis N, Rosen B, Dale AM. ''Automatically Parcellating the Human Cerebral Cortex.'' Cerebral Cortex (2004)} 14: 11-22. DOI: 10.1093/cercor/bbg087',...
    'Destrieux C, Fischl B, Dale A, Halgren E. ''Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature.'' Neuroimage (2010) 53(1):1-15.'});

lbls = readtable('src/destrieux_labels.csv');
labels = [cellfun(@(x1)(sprintf('L_%s', x1)), lbls.Var2(2:end), 'UniformOutput', false);...
    cellfun(@(x1)(sprintf('R_%s', x1)), lbls.Var2(2:end), 'UniformOutput', false)];
label_descriptions = [cellfun(@(x1)(sprintf('%s (left)', x1)), lbls.Var7(2:end), 'UniformOutput', false);...
    cellfun(@(x1)(sprintf('%s (right)', x1)), lbls.Var7(2:end), 'UniformOutput', false)];

pmap = {};
for study = {'bmrk5', 'paingen', 'spacetop'}
    for hemi = {'lh' 'rh'}
        pdata = fmri_data(which(sprintf('%s_%s_destrieux_%s.nii.gz', hemi{1}, study{1}, space_description)));
        
        pmap{end+1} = zeros(size(pdata.dat,1), height(lbls)-1);
        for i = 1:height(lbls)-1
            pmap{end}(:,i) = mean(pdata.dat == lbls.Var1(i+1),2);
        end 
    end
    % combine left and right hemispheres
    pmap{end-1} = cat(2,pmap{end-1:end});
    pmap(end) = [];
end

pdata = pdata.get_wh_image(1);

meanPmap = mean(cat(3,pmap{:}),3);
pdata.dat = meanPmap;

% renorm (mainly midline regions with lh/rh overlap)
total_p = sum(pdata.dat,2);
s = total_p(total_p > 1);
pdata.dat(total_p > 1,:) = pdata.dat(total_p>1,:)./s;

atlas_obj = atlas(pdata, ...
    'atlas_name', atlas_name, ...
    'labels', labels',...
    'label_descriptions', label_descriptions,...
    'space_description', space_description, ...
    'references', references, 'noverbose');
atlas_obj.probability_maps = sparse(atlas_obj.probability_maps);


% Check display
% -----------------------------------------------------------------------

dosave = true;

% Display with unique colors for each region:
atlas_obj.orthviews('unique',which([space_description '_T1_1mm.nii.gz']));


% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

cmap = [lbls.Var3, lbls.Var4, lbls.Var5];
cmap = cmap(2:end,:);
cmap = [cmap; cmap];
cmap = cmap(find(any(pdata.dat)),:)/255;

if dosave
   
    o2 = canlab_results_fmridisplay([], 'full hcp', ...
        'overlay', which(sprintf('%s_T1_1mm.nii.gz',space_description)));
    brighten(.6)
    
    o2 = montage(r, o2, 'indexmap', cmap);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

    
end

%% save object

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save(savename, 'atlas_obj');
    
end

%% write - this writes only the label image

if dosave
    
    savename = sprintf('%s_atlas_regions.img', atlas_name);
    atlas_obj.fullpath = fullfile(pwd, savename);
    write(atlas_obj,'overwrite');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names
labels = strrep(strrep(atlas_obj.labels,'&','_and_'),'-','_');
for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    arrayfun(@(x1)(set(x1,'FaceAlpha', .5)), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end


%% save nifti version
pdata.fullpath = sprintf('destrieux_%s.nii', space_description);
pdata.write('overwrite')
gzip(pdata.fullpath)
delete(pdata.fullpath)
writetable(table(atlas_obj.labels', atlas_obj.label_descriptions, 'VariableNames', {'labels', 'label_descriptions'}), 'destrieux_labels.csv');