function [match_result, matm] = match(tilde_p, tilde_q, Cx, Cy, ntrials, jsd_threshold, numCols, n_shuffles)

ucp = unique(Cx); ucq = unique(Cy);
nrowcluster1 = length(ucp);
nrowcluster2 = length(ucq);

separated_matrices_p = separate_clusters(Cx, tilde_p, numCols);
separated_matrices_q = separate_clusters(Cy, tilde_q, numCols);

matm = cell(n_shuffles, 1);
total_jsds = zeros(n_shuffles, 1);
valid_pair_counts = zeros(n_shuffles, 1);

for shuffle_idx = 1:n_shuffles
    fprintf('\n Attempting matching trial #%d:\n', shuffle_idx);

    order = randperm(nrowcluster2);
    matches = NaN(nrowcluster2, 3);
    used_sources = false(nrowcluster1, 1);

    for j = 1:nrowcluster2
        % If all source clusters have been used, no further matching is possible
        if all(used_sources)
            matches(j,1:2) = [NaN, ucq(order(j))];
            matches(j,3) = NaN;
            continue;
        end

        target_cluster = separated_matrices_q{order(j)};
        dis_values = Inf(nrowcluster1, 1);

        for i = 1:nrowcluster1
            if used_sources(i)
                continue;
            end
            source_cluster = separated_matrices_p{i};
            dis_values(i) = similarity(source_cluster, target_cluster, ntrials);
        end

        [min_jsd, best_match_index] = min(dis_values);

        if min_jsd < jsd_threshold
            matches(j,1) = ucp(best_match_index);         % matched source cluster
            matches(j,2) = ucq(order(j));                 % target cluster
            matches(j,3) = min_jsd;
            used_sources(best_match_index) = true;
        else
            matches(j,1:2) = [NaN, ucq(order(j))];         % no match, but JSD still meaningful
            matches(j,3) = min_jsd;
        end
    end

    matm{shuffle_idx} = matches;

    % Only include valid (non-NaN) JSD values in the total score
    valid_jsds = matches(~isnan(matches(:,3)), 3);
    total_jsds(shuffle_idx) = sum(valid_jsds);
    valid_pair_counts(shuffle_idx) = sum(~isnan(matches(:,1)));% Count valid matches based on source cluster column


    fprintf('Number of valid matched pairs: %d / %d\n', valid_pair_counts(shuffle_idx), nrowcluster2);
    fprintf('Total JSD distance: %.4f\n', total_jsds(shuffle_idx));
    disp(matches);
end

[~, best_idx] = min(total_jsds);
match_result = matm{best_idx};

fprintf('\n Final selected matching result is from trial #%d (Total JSD distance = %.4f):\n', best_idx, total_jsds(best_idx));
disp(match_result);

end
