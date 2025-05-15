function separated_matrices = separate_clusters(Cy, tilde_p, numCols)
% We only need to use the RNA part for matching
    tilde_p_subset = tilde_p(:, 1:numCols);
    unique_label = unique(Cy);
    separated_matrices = cell(length(unique_label), 1);
    
    % Iterate over each unique cluster label
    for i = 1:length(unique_label)
        % Create logical index
        cluster_mask = (Cy == unique_label(i));
        % Use logical index to extract corresponding rows
        separated_matrices{i} = tilde_p_subset(cluster_mask, :);
    end
end
