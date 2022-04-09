function obj = image_object_weighted_average(obj, indx1, indx2)

obj = replace_empty(obj);
obj.dat(obj.dat == 0) = NaN; % exclude NaNs, use existing valid data only where available

n1 = obj.metadata_table.N(indx1);
n2 = obj.metadata_table.N(indx2);

imgdat1 = obj.dat(:, indx1);
imgdat2 = obj.dat(:, indx2);

%imgdat = (imgdat1 .* n1 + imgdat2 .* n2) ./ (n1 + n2);

imgdat = (nansum([imgdat1 .* n1 imgdat2 .* n2]')') ./ (n1 + n2);

obj.dat(:, indx1) = imgdat;
obj.dat(:, indx2) = [];

% meta-data
k = size(obj.metadata_table, 2);
t = obj.metadata_table(indx1, :);
t2 = obj.metadata_table(indx2, :);

for i = 1:k
    
    if iscell(t{1, i})
        
    val1 = t{1, i};
    val2 = t2{1, i};
    t(1, i) = {['Average: ' val1{1} '_' val2{1}]};  end

end

t.n_male = t.n_male + t2.n_male;
t.n_female = t.n_female + t2.n_female;
t.N = t.N + t2.N;
t.age_years = mean([t.age_years t2.age_years]);

t.sd_age = mean([t.sd_age t2.sd_age]);
t.scan_length_min = mean([t.scan_length_min t2.scan_length_min]);

obj.metadata_table(indx1, :) = t;
obj.metadata_table(indx2, :) = [];

obj.image_names(indx1, :) = strrep(obj.image_names(indx1, :), '.nii', '.AVG');
obj.image_names(indx2, :) = [];

obj.fullpath(indx1, :) = strrep(obj.fullpath(indx1, :), '.nii', '.AVG');
obj.fullpath(indx2, :) = [];

obj.removed_images(indx2) = [];
obj.files_exist(indx2) = [];

obj.history{end+1} = sprintf('Averaged images %3.0f and %3.0f', indx1, indx2);

end

