function [tilde_p, cluster_p] = updateTildep_plus(p, Cx, Cz)

% Prevent empty cluster indices q
Cx = double(Cx(:));
Cz = double(Cz(:));

% Compute cluster_p: joint clustering distribution p(tilde_X, tilde_Z)
[i, j] = meshgrid(Cx, Cz);
subs = [i(:) j(:)];
val = p';
cluster_p = accumarray(subs, val(:), [max(Cx), max(Cz)], @sum, 0);
cluster_p = cluster_p / (sum(cluster_p(:)) + eps);

% Initialization
[nrow, ncol] = size(p);
tilde_p = zeros(nrow, ncol);

% Precompute components
sum_p_row = sum(p, 2);      % Total expression for each cell
sum_p_col = sum(p, 1);      % Total expression for each gene
sum_cp_row = sum(cluster_p, 2) + eps;   % Sum of clustering probabilities for each cell cluster
sum_cp_col = sum(cluster_p, 1)' + eps;  % Sum of clustering probabilities for each gene cluster

% Compute tilde_p in parallel
parfor i = 1:nrow
    r = Cx(i);
    row_val = sum_p_row(i);
    for j = 1:ncol
        c = Cz(j);
        tilde_p(i,j) = cluster_p(r,c) * ...
            row_val / sum_cp_row(r) * ...
            sum_p_col(j) / sum_cp_col(c);
    end
end

end
