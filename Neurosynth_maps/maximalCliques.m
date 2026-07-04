function cliques = maximalCliques(A)
% MAXIMALCLIQUES Find all maximal cliques in an undirected graph
% using the Bron–Kerbosch algorithm with pivoting.
%
% INPUT
%   A : NxN logical adjacency matrix (symmetric, no self-connections)
%
% OUTPUT
%   cliques : cell array of vectors, each containing node indices of a clique
%
% Example:
%   A = rand(10) > 0.7;
%   A = triu(A,1); A = A + A';
%   cliques = maximalCliques(A);

if ~islogical(A)
    A = A ~= 0;
end

% Remove self connections
A(1:size(A,1)+1:end) = 0;

n = size(A,1);
cliques = {};

% Initial sets
R = [];
P = 1:n;
X = [];

bronKerboschPivot(R, P, X);

% Nested function
    function bronKerboschPivot(R, P, X)
        if isempty(P) && isempty(X)
            cliques{end+1} = R; %#ok<AGROW>
            return;
        end

        % Choose pivot (heuristic: node with max neighbors in P ∪ X)
        PX = [P, X];
        if ~isempty(PX)
            degrees = sum(A(PX, P), 2);
            [~, idx] = max(degrees);
            u = PX(idx);
        else
            u = [];
        end

        % Candidates: P \ neighbors(u)
        if isempty(u)
            candidates = P;
        else
            candidates = setdiff(P, find(A(u,:)));
        end

        for v = candidates
            Nv = find(A(v,:));

            bronKerboschPivot([R, v], ...
                intersect(P, Nv), ...
                intersect(X, Nv));

            P = setdiff(P, v);
            X = [X, v];
        end
    end
end