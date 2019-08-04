input_var = ST_cleaned.cPDM;
mainanalysislabel = 'cPDM';

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

printhdr(' ');
printhdr(mainanalysislabel);
printhdr(' ');

%%
% Goal: Build an anonymous function called remove_nuisance that regresses out
% brain nuisance components related to brain measures of interest, but not
% pain/temperature:

% Correlations between raw extracted ROI values and potential nuisance components
% Question is: can we identify nuisance variables that are correlated with
% regional brain activity but not pain/temp?  And thus controlling for
% these, the relationship with pain becomes stronger. 
% We should consider whether we are dealing with raw or z-scored-within
% person/cleaned data.

% Create anonymous function for regressing an outcome on  multiple brain variables:
% (This builds in the particular brain variables in input_var as predictors)
% -------------------------------------------------------------------------

X = scale(input_var, 1);
pX = pinv(X);

get_brain_model_fits = @(y) X * pX * y;  % no intercept... ?? assuming all vars are z-scored...?
get_brain_model_varexplained = @(y) 100 .* (var(get_brain_model_fits(y)) ./ var(y));

%% Predict pain using non-residualized brain
% -------------------------------------------------------------------------
varexplained = get_brain_model_varexplained([ST_cleaned.rel_temp ST_cleaned.pain_rating]);

fprintf('%s variables explain %3.2f%% variance in temperature (relative values) and %3.2f%% variance in pain.\n', mainanalysislabel, varexplained(1), varexplained(2));


%% Create anonymous function for regressing an outcome on pain and temperature:
% (This builds in pain and temp as predictors)
% -------------------------------------------------------------------------

% behavior to nuisance 
X = scale([ST_cleaned.rel_temp ST_cleaned.pain_rating], 1);
pX = pinv(X);

paintempfit = @(y) X * pX * y;
paintempvarexplained = @(y) 100 .* (var(paintempfit(y)) ./ var(y));

%% Find nuisance variables (gray/white/csf) that relate to brain, but not pain/temp
% -------------------------------------------------------------------------
% we want to find nuisance covariates where v1/v2 is high

for i = 1:3
    
    v1 = get_brain_model_varexplained(ST_cleaned.gray_white_csf(:, i)); % brain to nuisance
    v2 = paintempvarexplained(ST_cleaned.gray_white_csf(:, i)); % behavior to nuisance
    fprintf('%s: brain %3.2f%% behavior %3.2f%% ratio (higher is better): %3.2f\n', ST_cleaned.graywhitecsf_labels{i}, v1, v2, v1/v2);
    
end

clear v1 v2

for i = 1:15
    
    v1(i) = get_brain_model_varexplained(ST_cleaned.gray_white_csf_components(:, i)); % brain to nuisance
    v2(i) = paintempvarexplained(ST_cleaned.gray_white_csf_components(:, i)); % behavior to nuisance
    fprintf('%s: brain %3.2f%% behavior %3.2f%% ratio (higher is better): %3.2f\n', ST_cleaned.graywhitecsfcomponents_labels{i}, v1(i), v2(i), v1(i)/v2(i));
    
end

wh = v1./v2 > 20;  % 20 times more predictive of brain than outcome. Good candidates for removal.

% Build an anonymous function called remove_nuisance that regresses out
% brain nuisance components related to brain measures of interest, but not
% pain/temperature:

X = scale(ST_cleaned.gray_white_csf_components(:, wh), 1);
X(:, end + 1) = 1; % intercept
remove_nuisance = @(y) y - X * pinv(X) * y;

%% Test variance explained in pain/temp after removing nuisance
% Fit models predicting pain with and without removal of CSF/etc components
% Multiple correlation: PDM1 without and with components removed

X = scale(input_var, 1);
X(:, end + 1) = 1;          % intercept

% Multiple regression table
[b, dev, stats] = glmfit(X(:, 1:end-1), ST_cleaned.pain_rating);
glm_table(stats, Xlabels);

disp(' ');

% Define anonymous functions
fit_brain = @(y) X * pinv(X) * y;                    
varexplained = @(y) var(fit_brain(y)) ./ var(y);

% Define anonymous functions with nuisance regression
X = remove_nuisance(X);
fit_brain_nonuisance = @(y) X * pinv(X) * y;                            % different X now...                 
varexplained_nonuisance = @(y) var(fit_brain_nonuisance(y)) ./ var(y);

var_explained_pain1 = varexplained(ST_cleaned.pain_rating);
var_explained_pain2 = varexplained_nonuisance(ST_cleaned.pain_rating);

fprintf('%s predicting pain: %3.2f%% variance explained without nuisance regression, %3.2f%% variance explained with; Ratio is %3.2f.\n', mainanalysislabel, 100*var_explained_pain1, 100*var_explained_pain2, var_explained_pain2./var_explained_pain1);

var_explained_temp1 = varexplained(ST_cleaned.rel_temp);
var_explained_temp2 = varexplained_nonuisance(ST_cleaned.rel_temp);

fprintf('%s predicting temperature: %3.2f%% variance explained without nuisance regression, %3.2f%% variance explained with; Ratio is %3.2f.\n', mainanalysislabel, 100*var_explained_temp1, 100*var_explained_temp2, var_explained_temp2./var_explained_temp1);

% cPDM predicting pain: 7.61% variance explained without nuisance regression, 7.64% variance explained with; Ratio is 1.00.
% cPDM predicting temperature: 6.20% variance explained without nuisance regression, 6.21% variance explained with; Ratio is 1.00.


