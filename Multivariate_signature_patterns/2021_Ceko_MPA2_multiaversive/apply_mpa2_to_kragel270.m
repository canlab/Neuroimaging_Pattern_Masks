[imgs, imgnames, allimgnames] = load_image_set('kragel18_alldata');

pexp = apply_multiaversive_mpa2_patterns(imgs, 'cosine_similarity');

pattnames = pexp.Properties.VariableNames;
pexpdat = table2array(pexp);

figure; plot(pexpdat); legend(pattnames);

set(gca, 'Xtick', 1:270, 'XTickLabel', allimgnames, 'XTickLabelRotation', 90);

studyid = imgs.dat_descrip.Studynumber;

studylist = studyid(1:15:end); % in non-sequential order that we want

%% Resort values into structure
% Fields are patterns from MPA2
% Cells are studies
% This format compatible with barplot_columns

PEXP = struct;
for i = 1:5, PEXP.(pattnames{i}) = []; end

for i = 1:length(studylist)
    
    wh = studyid == studylist(i);

    for j = 1:5
        % For each pattern
        
        PEXP.(pattnames{j}){i} = pexpdat(wh, j);
        
    end
    
end

%%
colors = [repmat({[1 0 0]}, 1, 2) repmat({[1 .5 0]}, 1, 2) repmat({[.7 .6 0]}, 1, 2) repmat({[0 0 1]}, 1, 2) repmat({[0 .5 1]}, 1, 2) repmat({[0 .6 .7]}, 1, 2) repmat({[0 1 0]}, 1, 2) repmat({[.5 1 0]}, 1, 2) repmat({[0 .7 .6]}, 1, 2)];

for i = 1:5
    
    create_figure(pattnames{i});
    
    barplot_columns(PEXP.(pattnames{i}), 'nofigure', 'noviolins', 'colors', colors);
    set(gca, 'XTickLabel', imgnames, 'XTickLabelRotation', 45);
    title(pattnames{i}, 'FontSize', 18);
    
    drawnow
    
end
