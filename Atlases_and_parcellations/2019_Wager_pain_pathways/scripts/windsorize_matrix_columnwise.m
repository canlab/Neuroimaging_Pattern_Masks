function [X, percent_adjusted] = windsorize_matrix_columnwise(X)
% [X, percent_adjusted] = windsorize_matrix_columnwise(X)

m = nanmean(X);
s = nanstd(X);
thr = 3 * s;        % standard deviations

[n, k] = size(X);

for i = 1:k
    
    wh = X(:, i) > m(i) + thr(i);
    
    percent_above(i) = 100 * sum(wh) ./ n;
    X(wh, i) = m(i) + thr(i);
    
    wh = X(:, i) < m(i) - thr(i);
    
    percent_below(i) = 100 * sum(wh) ./ n;
    X(wh, i) = m(i) - thr(i);
    
end

percent_adjusted = mean(percent_above + percent_below);

end

