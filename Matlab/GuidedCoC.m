function [Cx_truth_new, Cy, Cz, cluster_p, cluster_q, cluster_q0, matm, match_result, obj] = ...
    GuidedCoC(p, q, q0, Cx_truth, nrowcluster2, ncolcluster, iter, ...
    beta, alpha, ntrials, jsd_threshold, n_shuffles, numCols)

if isempty(gcp('nocreate')), parpool; end
rng(1);

Cx = mapLabels(Cx_truth);
Cx_truth_new = mapLabels(Cx_truth);
Cy = randsample(nrowcluster2, size(q, 1), true);
Cz = randsample(ncolcluster, size(p, 2), true);

[tilde_p, cluster_p]   = updateTildep_plus(p, Cx, Cz);
[tilde_q, cluster_q]   = updateTildep_plus(q, Cy, Cz);
[tilde_q0, cluster_q0] = updateTildep_plus(q0, Cy, Cz);

obj = zeros(iter, 1);



% Main loop
for t = 1:iter
    fprintf('Iteration %d\n', t);

    % Update Cy

    [Cy, dist0] = updateRowClustering_Y(q, tilde_q, q0, tilde_q0, Cy, alpha);
    [tilde_q, cluster_q]   = updateTildep_plus(q, Cy, Cz);
    [tilde_q0, cluster_q0] = updateTildep_plus(q0, Cy, Cz);
    % Update Cz

    [Cz, dist]  = updateColClustering_Z(p, q, q0, tilde_p, tilde_q, tilde_q0, Cz, beta, alpha);

    % Update tilde_p/q/q0

    [tilde_p, cluster_p]   = updateTildep_plus(p, Cx, Cz);
    [tilde_q, cluster_q]   = updateTildep_plus(q, Cy, Cz);
    [tilde_q0, cluster_q0] = updateTildep_plus(q0, Cy, Cz);

    % Save current loss and timestamp
    obj(t) = dist + dist0;
    
end


%  Matching  phase
fprintf(' Matching source and target cell clusters...\n');
[match_result, matm] = match(tilde_p, tilde_q, Cx, Cy, ntrials, jsd_threshold, numCols, n_shuffles);
disp(match_result);

end
