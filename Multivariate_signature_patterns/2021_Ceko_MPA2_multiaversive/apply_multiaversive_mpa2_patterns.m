% Given an fmri_data object with data, return pattern responses for 5 models in MPA2 multiaversive
% :Usage:
% ::
%
%     pexp_tabl = apply_multiaversive_mpa2_patterns(data_obj, [similarity_metric])
%
%     - Loads 5 MPA2 models and applies the linear patterns to data_obj
%     - Allows choice of similarity metric
%     - Adds intercept for each model if dot product (default) metric is
%     used
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Phil Kragel and Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **data_obj:**
%        An fmri_data object with one or more data images
%
% :Optional Inputs:
%
%   **similarity_metric:**
%        'dotproduct' (default). Will add intercepts 
%        'cosine_similarity'   . Will not add intercepts; they are no longer meaningful 
%
% :Outputs:
%
%   **pexp_tabl:**
%        A table of [images x 5] pattern expression values, with column
%        names indicating the name of the model.
%        Models are: {'general' 'mechanical' 'thermal' 'sound' 'visual'}
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :References:
%   Ceko et al. 2022 Nat Neurosci. Marta Ceko and Tor Wager.
%
% :See also:
%   - apply_all_signatures.m, apply_nps.m, apply_siips.m
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..
% Current version: 1.1
%
% Oct 2020:     Version 1.0 created by Phil Kragel
% Nov 4, 2020:  Version 1.1 created by Tor; add documentation, metric
%               options, intercepts.


% Version 
function pexp_tabl = apply_multiaversive_mpa2_patterns(dat, varargin)

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

similarity_metric = 'dotproduct';
intcpts = [0.1831233   0.0175146   0.1473605   0.0037268   0.0145212];

allowable_inputs = {'similarity_metric'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'cosine_similarity'
                
                similarity_metric = 'cosine_similarity';
                
            case 'correlation'
                
                similarity_metric = 'correlation';
                
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Zero out intercepts if they are not meaninful for chosen similarity metric
if strcmp(similarity_metric, 'cosine_similarity') || strcmp(similarity_metric, 'correlation')
    intcpts = [0 0 0 0 0];
end

% -------------------------------------------------------------------------
% LOAD PATTERNS (MODELS) AND APPLY THEM TO DATA
% -------------------------------------------------------------------------

    pats = load_image_set('mpa2');
    
    pexps = apply_mask(dat, pats, 'pattern_expression', similarity_metric);
    
    % Add intercepts for each model (or zeros for cosine_sim/corr)
    
    for i = 1:5
        pexps(:, i) = pexps(:, i) + intcpts(:, i);
    end
    
    pexp_tabl = array2table(pexps, 'VariableNames', {'general' 'mechanical' 'thermal' 'sound' 'visual'});
    

end
    
    

