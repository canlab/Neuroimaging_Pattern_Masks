function [SIG, sigtable] = apply_all_signatures(DATA_OBJ, varargin)
%
% Applies a series of 'signature' patterns to a set of image objects
% - Choice of image scaling method
% - Choice of similarity metric
% - Output is table (Matlab object type) format data in structure 
%
% :Usage:
% ::
%
%     SIG = apply_all_signatures(DATA_OBJ, ['similarity_metric', sim_metric, 'image_scaling', imgscaling, 'conditionnames', conditionnames, 'image_set', image_set_name)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2017 Tor Wager - 
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
%   **DATA_OBJ:**
%        Cell array of fmri_data objects, or single fmri_data object
%        Signatures are applied to each image in each object
%        Cells index 'conditions', which you can name
%
% :Optional Inputs:
%   **similarity_metric:**
%        Followed by similarity metric to use.  
%        Options are 'dotproduct' (default) 'cosine_similarity' 'correlation'
%        These options are passed into apply_mask, which uses them.
%
%   **image_scaling:**
%        Followed by image scaling metric to use.
%        Common options are 'none' (default), 'centerimages', 'zscoreimages', 'l2norm_images'
%        These options are passed into fmri_data.rescale, which uses them.
%
%   **conditionnames:**
%        Followed by cell array of names for each condition
%        WARNING: special characters will break this function, as these
%        names are used as variable names to define output tables.
%
%	%%:image_set:**
%       Followed by text indicating signature set to examine
%       See help load_image_set.m for options
%       e.g., 'npsplus'
%
%
% :Outputs:
%
%   **SIG:**
%        - Data structure with one field per named signature
%        - Each field contains a data table whose variables are named with
%        your conditionnames 
%        - Tables contain pattern response for images x conditions
%
%   **sigtable:**
%        - Cell vector of tables, one table per cell of DATA_OBJ
%        - Columns of each table are signatures
%
% :Examples:
% ::
%
%    [SIG, sigtable] = apply_all_signatures(DATA_OBJ, 'conditionnames', DAT.conditions);
%
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%    Created by Tor Wager, March 2017
%
%    Adapted for other image sets by pkragel, 5//16/2017
% ..

% -------------------------------------------------------------------------
% ..
%    DEFAULTS AND INPUTS
% ..
% -------------------------------------------------------------------------


similarity_metric = 'dotproduct';   % cosine_similarity correlation
image_scaling = 'none';             % 'l2norm_images' 'zscoreimages' 'centerimages'
image_set = 'npsplus';

% Make data into cell if not already
if ~iscell(DATA_OBJ), DATA_OBJ = {DATA_OBJ}; end

% number of conditions
k = length(DATA_OBJ);  

% Default names
for i = 1:k
    conditionnames{i} = sprintf('C%3.0f', i);
end

% Table format - tables for each condition
sigtable = cell(1, k);
[sigtable{:}] = deal(table());

% -------------------------------------------------------------------------
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'conditions', 'condition_names', 'object_names', 'conditionnames'}, conditionnames = varargin{i+1}; varargin{i+1} = [];
            case 'similarity_metric', similarity_metric = varargin{i+1}; varargin{i+1} = [];
            case 'image_scaling', image_scaling = varargin{i+1}; varargin{i+1} = [];
            case 'image_set', image_set= varargin{i+1}; varargin{i+1} = [];   
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% Load signatures - requires privately shared masks
% -------------------------------------------------------------------------
[signature_obj, signaturenames, imgnames] = load_image_set(image_set);

% resample_space: speed up later by avoiding reslicing
signature_obj = resample_space(signature_obj, DATA_OBJ{1});

% -------------------------------------------------------------------------
% Enforce formatting of names
% -------------------------------------------------------------------------

signaturenames = strrep(signaturenames, '-', '_');
signaturenames = strrep(signaturenames, ' ', '_');
signaturenames = strrep(signaturenames, '.', '');
signaturenames = strrep(signaturenames, '(', '');
signaturenames = strrep(signaturenames, ')', '');
signaturenames = strrep(signaturenames, '^', '_');
%mysigname = strrep(mysigname, ' ', '_');

% Remove blanks, common special chars, which will return errors
for i = 1:length(conditionnames), conditionnames{i} = strrep(conditionnames{i}, ' ', '_'); end
for i = 1:length(conditionnames), conditionnames{i} = strrep(conditionnames{i}, ',', '_'); end
for i = 1:length(conditionnames), conditionnames{i} = strrep(conditionnames{i}, '-', '_'); end

% -------------------------------------------------------------------------
% Initialize output structure
% -------------------------------------------------------------------------

SIG = struct('similarity_metric', similarity_metric, 'image_scaling', image_scaling);
SIG.signaturenames = signaturenames;
SIG.conditionnames = conditionnames;

% Get sizes for padding purposes, if different numbers of images in objects
for i = 1:k
    nimgs(i, 1) = size(DATA_OBJ{i}.dat, 2);
end
maxnimgvec = ones(max(nimgs), 1);

% -------------------------------------------------------------------------
% Normalize/scale objects, if requested
% -------------------------------------------------------------------------
switch image_scaling
    
    case 'none'
        % do nothing
        
    otherwise
        % use image_scaling as input to rescale method, for each object
        
        for i = 1:k
            
            DATA_OBJ{i} = rescale(DATA_OBJ{i}, image_scaling);
            
        end
        
end



% Main loop
% -------------------------------------------------------------------------
for j = 1:length(signaturenames)
    
    fprintf(1, '%s ', signaturenames{j});
    
    mysignature = get_wh_image(signature_obj, j);

    conditiontable = table();
    
    for i = 1:k
        
        myresponse = apply_mask(DATA_OBJ{i}, mysignature, 'pattern_expression', 'ignore_missing', similarity_metric);
        
        myresponse = padwithnan(myresponse, maxnimgvec, 1);
        
        conditiontable.(conditionnames{i}) = myresponse;

        mysigname = signaturenames{j};
        
        % tables for each condition
        varname = [mysigname '_' similarity_metric '_' image_scaling];

        sigtable{i}.(varname) = myresponse;
        
        sigtable{i}.Properties.VariableDescriptions{j} = sprintf('%s signature response. Similarity = %s, Data scaling = %s', mysigname, similarity_metric, image_scaling);
        
    end
    
    SIG.(mysigname) = conditiontable;
    

end

fprintf(1, '\n');
    
end % function
    
   