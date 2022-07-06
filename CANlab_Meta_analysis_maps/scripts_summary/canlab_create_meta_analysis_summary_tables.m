
all_n_tables = struct();

%%
cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2012_Denny_SOMA_self_other_meta')

dbname = 'SocialSelf';

all_n_tables.(dbname) = process_DB;

%%
cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2011_Agency_Meta_analysis')

dbname = 'Agency';

all_n_tables.(dbname) = process_DB;

%%

cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2014_BuhleSilvers_Reappraisal_Meta')

dbname = 'EmotionRegulation';

all_n_tables.(dbname) = process_DB;

%%

cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2018_Kraynak_Gianaros_immunemeta')
load('Kraynak_2018_Meta_DB_foci')
save SETUP DB
dbname = 'Neuroimmune';

all_n_tables.(dbname) = process_DB;

%%

cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2017_Ashar_Placebo_Review/Placebo_increases')

dbname = 'Placebo';

all_n_tables.(dbname) = process_DB;

%%

cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2015_Satpute_Barett_Emotion_Valence')
load Emotion_Meta_DB_all_3_30_13 DB
save SETUP DB

dbname = 'Emotion';

all_n_tables.(dbname) = process_DB;

%%

cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2015_Wager_Kang_etal_Emotion_Meta_BSPP')
load Emotion_Meta_DB_5Emotions_3_30_13 DB
save SETUP DB

dbname = 'Emotion_Categories';

all_n_tables.(dbname) = process_DB;

%%

cd('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/scripts_summary')

save meta_analysis_summary_tables all_n_tables


%% WORK IN PROGRESS: see what we can recover from old databases
% Triage this for now unless someone really wants it
% dbs = load('/Users/f003vz1/Documents/GitHub/Canlab_MKDA_MetaAnalysis/meta_analysis_master_file.mat')
% dbs.INHDB
% DB = Meta_Setup(dbs.INHDB);


%%


function N_by_year_table = process_DB

load SETUP DB
DB = Meta_tables_of_variables(DB);
save SETUP -append DB

N_by_year_table = varfun(@mean, DB.TABLES.study_table, 'GroupingVariables','Year', 'InputVariables', 'mean_Subjects', 'OutputFormat','table');
ste_table = varfun(@ste, DB.TABLES.study_table, 'GroupingVariables','Year', 'InputVariables', 'mean_Subjects', 'OutputFormat','table');

N_by_year_table = join(N_by_year_table, ste_table);

create_figure('Sample size'); 
plot(DB.TABLES.study_table.Year, DB.TABLES.study_table.mean_Subjects, 'ko');
plot(N_by_year_table.Year, N_by_year_table.mean_mean_Subjects, 'LineWidth', 3, 'Color', 'b');

xlabel('Year')
ylabel('Mean sample size per contrast')

end
