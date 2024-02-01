function [siips_values, image_names, data_objects, siipspos_exp_by_region, siipsneg_exp_by_region, clpos, clneg] = apply_siips(input_images, varargin)
%
% Applies the SIIPS1 signuature pattern to a set of brain images (nii or img)
% - requires SIIPS1 mask and related extensions on path with standard name(s)
% - takes input in multiple forms, including wildcard, filename list, and cell of fmri_data objects
% - this can also be done easily using apply_mask.m for any mask/signature,
%   but apply_siips makes it easier.
% - apply_siips also extracts local pattern responses from SIIPS subregions
%
% Usage:
% -------------------------------------------------------------------------
% [outputs] = function_name(list inputs here, [optional inputs])
%
% [siips_values, image_names, data_objects, siipspos_exp_by_region, siipsneg_exp_by_region, siipspos, siipsneg] = apply_siips(input_images)
%
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2014 Tor Wager
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
%
% Inputs:
% -------------------------------------------------------------------------
% input_images           multiple formats:
%                        - fmri_data object (or other image_vector object)
%                        - fmri_data objects (or other image_vector objects) in cell { }
%                        - list of image names
%                        - cell { } with multiple lists of image names
%                        - image wildcard (NOTE: uses filenames.m, needs
%                        system compatibility)
%                        - cell { } with multiple image wildcards
% Optional inputs:
% 'noverbose'           suppress screen output
%                       not recommended for first run, as verbose output
%                       prints info about missing voxels in each image
% 'notables'            suppress text table output with images and values
%
%                       Note:Similarity metric passed through to canlab_pattern_similarity
% 'cosine_similarity'   Use cosine similarity measure instead of dot product
% 'correlation'         Use correlation measure instead of dot product
% 'donorm'              L1 norm the pattern expression for all regions

% Outputs:
% -------------------------------------------------------------------------
% siips_values            cell { } with siips values for each image
% image_names
% data_objects
% siipspos_exp_by_region  local pattern expression in subregions defined by siips pos mask
% siipsneg_exp_by_region  local pattern expression in subregions defined by siips neg mask
%                       These are matrices of observations x regions,
%                       placed in a cell array with one cell per image
%                       list/fmri_data object entered as input
% clpos                siips positive subregions, region object. Names of regions in cl(x).shorttitle
% clneg                siips negative subregions, region object. Names of regions in cl(x).shorttitle
%
% Examples:
% -------------------------------------------------------------------------
% % Enter a pre-specified list of images:
% input_images = filenames('image_data/Pain_Sub*ANP_001.img', 'char', 'absolute');
% [siips_values, image_names, data_objects] = apply_siips(input_images, 'noverbose');
%
% % % Enter wildcard (requires system compatibility):
% [siips_values, image_names, data_objects] = apply_siips('image_data/Pain_Sub*ANP_001.img', 'noverbose');
%
% % Enter pre-defined objects:
% data_object = load_image_set('emotionreg');
% [siips_values, image_names] = apply_siips(data_object, 'noverbose');
%
% % Enter series of wildcards:
% wcards = {'Pain_Sub*AP_001.img' ...
%           'Pain_Sub*ANP_001.img' ...
%           'Pain_Sub*EP_001.img' ...
%           'Pain_Sub*ENP_001.img' ...
%           'Pain_Sub*AP_002.img' ...
%           'Pain_Sub*ANP_002.img' ...
%           'Pain_Sub*EP_002.img' ...
%           'Pain_Sub*ENP_002.img' ...
%           'Pain_Sub*AP_003.img' ...
%           'Pain_Sub*ANP_003.img' ...
%           'Pain_Sub*EP_003.img' ...
%           'Pain_Sub*ENP_003.img'};
%
% [siips_values, image_names, data_objects] = apply_siips(wcards);
%
% Extract and return individual local pattern expression values:
% [siips_values, ~, ~, siipspos_exp_by_region, siipsneg_exp_by_region] = apply_siips(data_objects);
%
% Example: Extract subregions and plot data in positive subregions with names:
% If you have the canlab_core and the new repository on your path with subfolders, 
% the code below provides a complete example of extracting the SIIPS subregions and 
% creating a bar plot showing pattern responses for a subset of the SIIPS subregions:
%
% data_object = load_image_set('emotionreg');
% [siips_values, image_names, data_objects, siipspos_exp_by_region, siipsneg_exp_by_region, clpos, clneg] = apply_siips(data_object);
% posnames = {clpos.shorttitle};
% barplot_columns(siipspos_exp_by_region{1}, 'names', posnames);
%
% See also:
% apply_mask.m

% Programmers' notes:
% Created by Tor Wager, Aug 2014

% INITIALIZE VARIABLES AND DEFAULTS
% -------------------------------------------------------------------------

verbose = 1;
verbosestr = 'verbose';
dotables = 1;
normstr='nonorm'; % default is No L1 Norm
mask = which('nonnoc_v11_4_137subjmap_weighted_mean.nii');

if nargout > 3
    % for individual subregions
    
    [siipspos, siipsneg, siips_pos_regions, siips_neg_regions, names_pos, names_neg, r] = load_siips_subregions();
    
end

if isempty(mask)
    error('Image ''weights_NSF_grouppred_cvpcr.img'' not found on path.');
end

[data_objects, image_names, data_objects] = deal({});

% SET UP OPTIONAL INPUTS
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'noverbose', verbose = 0; verbosestr = 'noverbose'; dotables = 0;
            case 'notables', dotables = 0;
            case 'donorm', normstr='donorm';
            case 'nonorm', normstr='nonorm';

            case {'cosine_similarity', 'correlation'}
                % do nothing here - these will be passed in to apply_mask
                % and applied in that function.
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% PARSE INPUT TYPE
% Return data_objects{} and image_names{}
% -------------------------------------------------------------------------

if isa(input_images, 'image_vector')
    % - fmri_data object (or other image_vector object)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    data_objects{1} = input_images;
    image_names{1} = data_objects{1}.image_names;
    clear input_images
    
elseif iscell(input_images) && isa(input_images{1}, 'image_vector')
    % - fmri_data objects (or other image_vector objects) in cell { }
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    data_objects = input_images;
    for i = 1:length(data_objects)
        image_names{i} = data_objects{i}.image_names;
    end
    clear input_images
    
elseif ischar(input_images) && ~any(input_images(:) == '*')
    % - list of image names
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    data_objects{1} = fmri_data(input_images, [], verbosestr);
    image_names{1} = input_images;
    clear input_images
    
elseif iscell(input_images) && ~any(input_images{1}(:) == '*')
    % - cell { } with multiple lists of image names
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for i = 1:length(input_images)
        data_objects{i} = fmri_data(input_images{i}, [], verbosestr);
        image_names{i} = input_images{i};
    end
    
elseif ischar(input_images) && size(input_images, 1) == 1 && any(input_images(:) == '*')
    % - image wildcard
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    image_names{1} = get_image_names(input_images);
    data_objects{1} = fmri_data(image_names{1}, [], verbosestr);
    clear input_images
    
elseif iscell(input_images) && size(input_images{1}, 1) == 1 && any(input_images{1}(:) == '*')
    % - cell { } with multiple image wildcards
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    for i = 1:length(input_images)
        image_names{i} = get_image_names(input_images{i});
        data_objects{i} = fmri_data(image_names{i}, [], verbosestr);
    end
    clear input_images
else
    error('Unknown operation mode or input type - check inputs and code.');
    
end % setup input type



% GET ALL siips VALUES
% -------------------------------------------------------------------------

% if nargout > 3
%     % for individual subregions
%     % Prep data. Sample data to siips mask so that number of regions does
%     % not change as a function of the input data. SEE BELOW.
%     siipspos = resample_space(fmri_data(siipspos, [], verbosestr), data_objects{1});
%     siipsneg = resample_space(fmri_data(siipsneg, [], verbosestr), data_objects{1});
% end

for i = 1:length(data_objects)
    
    % Pass in cosine similarity or correlation keywords
    if verbose
        siips_values{i} = apply_mask(data_objects{i}, mask, 'pattern_expression', varargin{:}); %, 'ignore_missing');
        
    else
        siips_values{i} = apply_mask(data_objects{i}, mask, 'pattern_expression', 'ignore_missing', varargin{:});
    end
    
    
    if nargout > 3
        
        % for individual subregions
        % Prep data. Sample data to siips mask so that number of regions does
        % not change as a function of the input data.
        data_objects{i} = resample_space(data_objects{i}, siipspos);
        
        % for individual subregions
        % no norm so that we preserve relative contributions to whole
        
        % Similarity metric passed in, passed through to canlab_pattern_similarity
        clpos = extract_roi_averages(data_objects{i}, siipspos, 'pattern_expression', 'contiguous_regions', normstr, verbosestr, varargin{:});
        clneg = extract_roi_averages(data_objects{i}, siipsneg, 'pattern_expression', 'contiguous_regions', normstr, verbosestr, varargin{:});

        % add names
        for j = 1:length(clpos)
            
            clpos(j).shorttitle = siips_pos_regions(j).shorttitle;
            
        end
        
        for j = 1:length(clneg)
            
            clneg(j).shorttitle = siips_neg_regions(j).shorttitle;
            
        end
        
        % save data from individual regions
        
        siipspos_exp_by_region{i} = cat(2, clpos.dat);
        siipsneg_exp_by_region{i} = cat(2, clneg.dat);
        
    end
    %         % re-sort data from individual regions
    %
    %     for j = 1:length(clpos)
    %         siipspos_exp_by_region{j}(:, i) = clpos(j).dat;
    %     end
    %
    %     for j = 1:length(clneg)
    %         siipsneg_exp_by_region{j}(:, i) = clneg(j).dat;
    %     end
    
end

% PRINT TABLE OF siips VALUES
% -------------------------------------------------------------------------
if verbose || dotables
    for i = 1:length(data_objects)
        
        if length(image_names) < i || isempty(image_names{i})
            image_names{i} = repmat('Unknown image names', size(siips_values{i}, 1), 1);
        elseif size(image_names{i}, 1) < size(siips_values{i}, 1)
            image_names{i} = repmat(image_names{i}(1, :), size(siips_values{i}, 1), 1);
        end
        
        Image = image_names{i};
        SIIPS1 = siips_values{i};
        mytable = table(Image, SIIPS1);
        
        if nargout > 3
            
            for j = 1:length(clpos)
                
                myname = format_strings_for_legend(clpos(j).shorttitle);
                myname = myname{1};
                myname = strrep(myname, ' ', '_');
                myname = strrep(myname, '(', '_');
                myname = strrep(myname, ')', '_');
                myname = strrep(myname, '/', '_');
                myname = strrep(myname, '\', '_');
                mytable.(myname) = siipspos_exp_by_region{i}(1:height(mytable), j);
                
            end
            
            for j = 1:length(clneg)
                
                myname = format_strings_for_legend(clneg(j).shorttitle);
                myname = myname{1};
                myname = strrep(myname, ' ', '_');
                myname = strrep(myname, '(', '_');
                myname = strrep(myname, ')', '_');
                myname = strrep(myname, '/', '_');
                myname = strrep(myname, '\', '_');
                mytable.(myname) = siipsneg_exp_by_region{i}(1:height(mytable), j);
                
            end
            
            disp(mytable)
            
        else
            
            print_matrix(siips_values{i}, {['siips values for series ' num2str(i)]}, cellstr(image_names{i}));
            %mytable = table(Image, SIIPS1)
            
        end
        
        
    end
end

end  % main function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function image_names = get_image_names(input_images)

wcard = input_images;
image_names = filenames(wcard, 'char', 'absolute');
if isempty(image_names)
    error(['No images found for wildcard ' wcard])
end

end % function