function [Cy, dist0] = updateRowClustering_Y(q, tilde_q, q0, tilde_q0, Cy, alpha)

% Initialize outputs
dist0 = NaN;

% compute q(Z|y) and q0(Z|y)
qZy  = q  ./ max(sum(q,2), eps);
q0Zy = q0 ./ max(sum(q0,2), eps);

% compute tilde_q(Z|tilde_y) and tilde_q0(Z|tilde_y)
tilde_qZy  = tilde_q  ./ max(sum(tilde_q,2), eps);
tilde_q0Zy = tilde_q0 ./ max(sum(tilde_q0,2), eps);

% Number of clusters (should be set to target number of clusters, not the number of unique values)

k = max(Cy);  
n = size(qZy, 1);  % 样本数

% Preallocate
tilde_qZty  = zeros(k, size(tilde_q,2));
tilde_q0Zty = zeros(k, size(tilde_q0,2));

% Adjust cluster labels to be 1~k to avoid out-of-bounds or empty clusters in accumarray

Cy = mapLabels(Cy);

% Use accumarray to aggregate and add eps to avoid clusters with all zeros

for i = 1:size(tilde_q,2)
    tmp = accumarray(Cy, tilde_qZy(:,i), [k,1], @sum, eps);
    tilde_qZty(:,i) = tmp;
end
for i = 1:size(tilde_q0,2)
    tmp = accumarray(Cy, tilde_q0Zy(:,i), [k,1], @sum, eps);
    tilde_q0Zty(:,i) = tmp;
end

% Initialize KL divergence matrix

temp = zeros(k, n);  

% Use parfor to compute KL divergence for each sample
parfor r = 1:n
    temp_r = zeros(k,1);
    for rc = 1:k
        temp_r(rc) = KLDiv(qZy(r,:), tilde_qZty(rc,:)) + ...
                     alpha * KLDiv(q0Zy(r,:), tilde_q0Zty(rc,:));
    end
    temp(:,r) = temp_r;
end

% Find the cluster with the minimum KL divergence
[dist_vec, best_assign] = min(temp);


if length(best_assign) ~= n
    error('Length mismatch: best_assign (%d) vs sample size (%d)', length(best_assign), n);
end

Cy = best_assign(:); 
dist0 = sum(dist_vec);


end
