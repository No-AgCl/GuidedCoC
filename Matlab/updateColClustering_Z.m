function [Cz, dist] = updateColClustering_Z(p, q, q0, tilde_p, tilde_q, tilde_q0, Cz, beta, alpha)

% Normalize Cz cluster labels 
Cz = mapLabels(Cz);
k = length(unique(Cz));  % Number of feature clusters

% Conditional distribution (Z|x) or (Z|y)
pXz  = p  ./ max(sum(p,1), eps);
qYz  = q  ./ max(sum(q,1), eps);
q0Yz = q0 ./ max(sum(q0,1), eps);

% Cluster-level distributions (Z|cluster z)
tilde_pXtz  = zeros(size(tilde_p,1), k);
tilde_qYtz  = zeros(size(tilde_q,1), k);
tilde_q0Ytz = zeros(size(tilde_q0,1), k);

% Accumulate soft assignments
for i = 1:size(tilde_p,1)
    tilde_pXtz(i,:) = accumarray(Cz, tilde_p(i,:)', [k,1], @sum, 0)';
end
for i = 1:size(tilde_q,1)
    tilde_qYtz(i,:) = accumarray(Cz, tilde_q(i,:)', [k,1], @sum, 0)';
end
for i = 1:size(tilde_q0,1)
    tilde_q0Ytz(i,:) = accumarray(Cz, tilde_q0(i,:)', [k,1], @sum, 0)';
end

% Prior for each feature
pz  = sum(p,1); 
qz  = sum(q,1); 
q0z = sum(q0,1);
num_z = size(pXz,2);  % number of features

% Initialize loss matrices
temp_X  = zeros(k, num_z);
temp_Y  = zeros(k, num_z);
temp_Y0 = zeros(k, num_z);

% Parallelized computation of KL divergences
parfor z = 1:num_z
    for zc = 1:k
        temp_X(zc, z)  = pz(z)  * KLDiv(pXz(:,z)'  , tilde_pXtz(:,zc)'  + eps);
        temp_Y(zc, z)  = qz(z)  * KLDiv(qYz(:,z)'  , tilde_qYtz(:,zc)'  + eps);
        temp_Y0(zc, z) = q0z(z) * KLDiv(q0Yz(:,z)' , tilde_q0Ytz(:,zc)' + eps);
    end
end

% Weighted objective
temp = temp_Y + beta * temp_X + alpha * temp_Y0;

% Update Cz: assign each feature to the best cluster
[mindist, Cz] = min(temp);
Cz = Cz';  % Ensure column vector


if length(Cz) ~= size(p, 2)
    error('Length of Cz (%d) does not match number of features (%d)', length(Cz), size(p, 2));
end

dist = sum(mindist);

clearvars -except Cz dist
end
