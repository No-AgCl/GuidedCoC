function dis = similarity(tilde_p, tilde_q, num_trials)
    % Get the number of rows for both matrices
    row_p = size(tilde_p, 1);
    row_q = size(tilde_q, 1);

    % Select the smaller number of rows
    row_num = min(row_p, row_q);

    % Use parfor to compute JSD across multiple trials in parallel

    jsd_values = zeros(num_trials, 1);
    parfor trial = 1:num_trials
        % Randomly select rows
        if row_p > row_num
            indices_p = randperm(row_p, row_num);
            new_tilde_p = tilde_p(indices_p, :);
        else
            new_tilde_p = tilde_p;
        end

        if row_q > row_num
            indices_q = randperm(row_q, row_num);
            new_tilde_q = tilde_q(indices_q, :);
        else
            new_tilde_q = tilde_q;
        end

        % Ensure one-to-one matching by row and compute JSD row by row

        jsd_single = zeros(row_num, 1);
        for i = 1:row_num
            jsd_single(i) = JSDiv(new_tilde_p(i,:), new_tilde_q(i,:));
        end

        jsd_values(trial) = mean(jsd_single);
    end

    % Return the average JSD over multiple trials
    dis = mean(jsd_values);
end
