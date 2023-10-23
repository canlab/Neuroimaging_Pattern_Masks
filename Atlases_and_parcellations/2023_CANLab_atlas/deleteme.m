% bstem_atlas_old = bstem_atlas
% biancia_old = biancia
bstem_atlas = bstem_atlas_old;

bstem_atlas = biancia.select_atlas_subset(biancia_regions).merge_atlases(bstem_atlas, 'noreplace');

target_atlas = atlas_obj;
uniq_roi = unique(target_atlas.remove_empty.dat);
target_atlas = target_atlas.replace_empty;
pmap = zeros(size(target_atlas.dat,1), length(uniq_roi));
for i = 1:length(uniq_roi)
    pmap(target_atlas.dat == i, i) = 0.35;
end
target_atlas.dat = pmap;
target_atlas.fullpath = '/tmp/bstem.nii';
target_atlas.write('overwrite')
