function [unique_cliques, out_obj, membership, stats] = assign_unique_cliques_from_maximal_cliques(obj, R, cliques, thr, names)
% ASSIGN_UNIQUE_CLIQUES_FROM_MAXIMAL_CLIQUES
%
% Start from a set of maximal cliques (e.g., output of maximalCliques.m),
% and assign each variable/image to exactly one clique.
%
% Rule:
%   For each variable i, consider each current clique G such that
%   corr(i,j) >= thr for all j in G\{i}. Among eligible cliques, assign i
%   to the clique with the highest minimum correlation to its members:
%
%       score(i,G) = min( R(i, G\{i}) )
%
%   If no current clique can validly accept i, assign i to its own
%   singleton clique.
%
% Iterate until convergence, pruning cliques that disappear.
%
% Then average maps in obj within each final unique clique, using
% get_wh_image(obj, idx), and store the averages in out_obj.
%
% Also creates:
%   out_obj.metadata_table
%
% with variables:
%   - original_indices : original image indices in each clique
%   - original_names   : concatenated names from input "names"
%   - clique_size      : size of each clique
%
% INPUTS
%   obj     : fmri_data-like object with one image per column
%   R       : NxN symmetric correlation/similarity matrix
%   cliques : cell array of index vectors, e.g. output of maximalCliques(A)
%   thr     : threshold, e.g. 0.6
%   names   : Nx1 or 1xN cell array / string array of original image names
%
% OUTPUTS
%   unique_cliques : cell array of final nonempty exclusive groups
%   out_obj        : object with one averaged image per group
%   membership     : Nx1 integer vector, membership(i) = group index
%   stats          : summary struct

n = size(R, 1);

if size(R, 2) ~= n
    error('R must be square.');
end

if ~iscell(cliques)
    error('cliques must be a cell array of index vectors.');
end

if isstring(names)
    names = cellstr(names(:));
elseif ischar(names)
    names = cellstr(names);
elseif ~iscell(names)
    error('names must be a cell array or string array.');
end

names = names(:);

if numel(names) ~= n
    error('names must have one entry per original image/variable.');
end

R = double(R);
R = (R + R') ./ 2;
R(1:n+1:end) = 1;

% Clean and standardize clique input
cliques = cliques(:);
for k = 1:numel(cliques)
    cliques{k} = unique(cliques{k}(:)'); %#ok<AGROW>
    cliques{k} = cliques{k}(cliques{k} >= 1 & cliques{k} <= n);
end
cliques = cliques(~cellfun(@isempty, cliques));

% If any variable is not present in any clique, seed a singleton clique for it
covered = false(n,1);
for k = 1:numel(cliques)
    covered(cliques{k}) = true;
end
missing_vars = find(~covered);
for ii = 1:numel(missing_vars)
    cliques{end+1} = missing_vars(ii); %#ok<AGROW>
end

% Initialize current cliques from input
unique_cliques = cliques;

max_iter = 200;
iter = 0;
converged = false;

while ~converged && iter < max_iter
    iter = iter + 1;

    old_unique_cliques = canonicalize_cliques(unique_cliques);

    % Current candidate groups
    groups = unique_cliques;
    ng = numel(groups);

    % Reassign every variable to exactly one current group
    membership = zeros(n,1);

    for i = 1:n
        best_group = [];
        best_score = -Inf;
        best_size  = -Inf;

        for g = 1:ng
            members = groups{g};
            others = members(members ~= i);

            if isempty(others)
                % singleton candidate: always eligible, weakest possible score
                eligible = true;
                score = -Inf;
                group_size = 1;
            else
                rvals = R(i, others);
                eligible = all(rvals >= thr);

                if eligible
                    score = min(rvals);    % maximin criterion
                    group_size = numel(members);
                end
            end

            if eligible
                % Prefer:
                % 1) larger minimum correlation
                % 2) larger group if tied
                % 3) first encountered if still tied
                if (score > best_score) || ...
                   (score == best_score && group_size > best_size)
                    best_group = g;
                    best_score = score;
                    best_size = group_size;
                end
            end
        end

        % If no current clique can accept i, create singleton
        if isempty(best_group)
            groups{end+1} = i; %#ok<AGROW>
            membership(i) = numel(groups);
        else
            membership(i) = best_group;
        end
    end

    % Rebuild exclusive cliques from membership
    new_cliques = cell(max(membership), 1);
    for g = 1:max(membership)
        new_cliques{g} = find(membership == g)';
    end

    % Prune empty/disappeared cliques
    new_cliques = new_cliques(~cellfun(@isempty, new_cliques));

    % Ensure each resulting group remains a clique under thr
    new_cliques = enforce_clique_constraint(new_cliques, R, thr);

    % Canonicalize for comparison
    unique_cliques = canonicalize_cliques(new_cliques);

    converged = isequal(unique_cliques, old_unique_cliques);
end

% Final membership vector
membership = zeros(n,1);
for g = 1:numel(unique_cliques)
    membership(unique_cliques{g}) = g;
end

% Build averaged output object
k = numel(unique_cliques);
out_obj = get_wh_image(obj, 1:k);

for g = 1:k
    idx = unique_cliques{g};
    this_obj = get_wh_image(obj, idx);
    out_obj.dat(:, g) = mean(this_obj.dat, 2, 'omitnan');
end

% Stats
stats = struct();
stats.iterations = iter;
stats.converged = converged;
stats.threshold = thr;
stats.n_items = n;
stats.n_groups = k;
stats.group_sizes = cellfun(@numel, unique_cliques(:));
stats.within_group_min_r = nan(k,1);
stats.within_group_mean_r = nan(k,1);

for g = 1:k
    idx = unique_cliques{g};
    if numel(idx) == 1
        stats.within_group_min_r(g) = 1;
        stats.within_group_mean_r(g) = 1;
    else
        subR = R(idx, idx);
        mask = tril(true(size(subR)), -1);
        vals = subR(mask);
        stats.within_group_min_r(g) = min(vals);
        stats.within_group_mean_r(g) = mean(vals, 'omitnan');
    end
end

% Optional labels if supported
label_strings = cell(k,1);
for g = 1:k
    label_strings{g} = sprintf('unique_clique_%03d_n%d', g, numel(unique_cliques{g}));
end

try
    out_obj.fullpath = label_strings;
end
try
    out_obj.image_names = label_strings;
end

% Build metadata_table
original_indices = cell(k,1);
original_names   = cell(k,1);
clique_size      = stats.group_sizes(:);

for g = 1:k
    idx = unique_cliques{g};
    original_indices{g} = idx(:)';
    original_names{g} = strjoin(names(idx), ' | ');
end

out_obj.metadata_table = table(original_indices, original_names, clique_size);

end


function cliques = enforce_clique_constraint(cliques, R, thr)
% Split any invalid group into a retained valid core plus singleton leftovers.
% Greedy but safe: final output groups all satisfy clique condition.

changed = true;

while changed
    changed = false;
    new_cliques = {};

    for c = 1:numel(cliques)
        idx = unique(cliques{c}(:)');

        if numel(idx) <= 1
            new_cliques{end+1} = idx; %#ok<AGROW>
            continue;
        end

        subR = R(idx, idx);
        offdiag_ok = subR >= thr | eye(numel(idx));

        if all(offdiag_ok(:))
            new_cliques{end+1} = idx; %#ok<AGROW>
        else
            % Remove the worst node: smallest minimum within-group r
            minr = inf(1, numel(idx));
            for j = 1:numel(idx)
                others = setdiff(1:numel(idx), j);
                minr(j) = min(subR(j, others));
            end
            [~, worst] = min(minr);

            kept = idx;
            removed = idx(worst);
            kept(worst) = [];

            if ~isempty(kept)
                new_cliques{end+1} = kept; %#ok<AGROW>
            end
            new_cliques{end+1} = removed; %#ok<AGROW>

            changed = true;
        end
    end

    cliques = new_cliques(~cellfun(@isempty, new_cliques));
end

% Merge duplicates if any arose
cliques = canonicalize_cliques(cliques);

end


function cliques = canonicalize_cliques(cliques)
% Sort members within each clique, remove empties, remove duplicates,
% and sort list stably for comparisons.

cliques = cliques(:);
cliques = cliques(~cellfun(@isempty, cliques));

for i = 1:numel(cliques)
    cliques{i} = unique(sort(cliques{i}(:)')); %#ok<AGROW>
end

% Remove duplicates
keys = cellfun(@(x) sprintf('%d_', x), cliques, 'UniformOutput', false);
[~, ia] = unique(keys, 'stable');
cliques = cliques(ia);

% Stable ordering
mins = cellfun(@min, cliques);
sizes = cellfun(@numel, cliques);
T = table((1:numel(cliques))', mins(:), sizes(:), 'VariableNames', {'idx','mins','sizes'});
T = sortrows(T, {'mins','sizes','idx'});
cliques = cliques(T.idx);

end