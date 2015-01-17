num_cuboids = length(joint_mean);

for i = 1:num_cuboids
    for j = 1:num_cuboids
        if i == j; continue; end;
		
        % joint_mean: Transposed.
        M = [joint_mean{i}{j}'; joint_inv_cov{i}{j}];
        csvwrite(['joint_normal_', num2str(i-1), '_',...
            num2str(j-1), '.csv'], M);
    end
end
