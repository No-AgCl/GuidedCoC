function CH_index = CH(X, num_row_clusters, num_col_clusters, mode)
    % 计算 CH 指数
    % mode = 'row' -> 计算最佳行聚类数（细胞聚类数）
    % mode = 'col' -> 计算最佳列聚类数（基因/Peak 聚类数）

    if num_clusters <= 0
        error('簇数量必须大于 0');
    end
    if num_clusters == 1
        CH_index = NaN;
        return;
    end
    
    % **选择 CoC 计算的模式**
    if strcmp(mode, 'row')
        % **计算最佳行聚类数（细胞聚类数）**
        [Cx, ~, ~, ~] = CoC(X, num_row_clusters, 10, 15);  % 10是默认列聚类数
        cluster_labels = Cx;  % 细胞的聚类标签
    elseif strcmp(mode, 'col')
        % **计算最佳列聚类数（基因/Peak 聚类数）**
        [~, Cz, ~, ~] = CoC(X, num_row_clusters, num_col_clusters, 15);  % 10是默认行聚类数
        cluster_labels = Cz;  % 特征（基因或 Peak）的聚类标签
    else
        error('mode 参数错误，应为 "row" 或 "col"');
    end

    % **计算 CH 指数**
    centers = zeros(num_clusters, size(X, 2));
    for i = 1:num_clusters
        centers(i, :) = mean(X(cluster_labels == i, :), 1);
    end
    
    SSB = 0;
    SSW = 0;
    total_mean = mean(X, 1);
    
    for i = 1:num_clusters
        cluster_data = X(cluster_labels == i, :);
        SSB = SSB + size(cluster_data, 1) * sum((mean(cluster_data, 1) - total_mean).^2);
        SSW = SSW + sum(pdist2(cluster_data, mean(cluster_data, 1)).^2);
    end
    
    CH_index = (SSB / (num_clusters - 1)) / (SSW / (size(X, 1) - num_clusters));
end



