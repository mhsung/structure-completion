num_cuboids = length(conditional_mean_A);

for i = 1:num_cuboids
    for j = 1:num_cuboids
        if i == j; continue; end;
        % conditional_mean_A, conditional_mean_b: Transposed.
        M = [conditional_mean_A{i}{j}'; conditional_mean_b{i}{j}';...
            conditional_inv_cov{i}{j}];
        csvwrite(['conditional_normal_', num2str(i-1), '_',...
            num2str(j-1), '.csv'], M);
    end
end
