num_cuboids = length(joint_mean);

for i = 1:num_cuboids
    for j = 1:num_cuboids
        if i == j; continue; end;
        
        fid = fopen(['joint_normal_', num2str(i-1), '_', num2str(j-1), '.dat'], 'w');
        
        fwrite(fid, size(joint_mean{i}{j}, 1), 'int16');
        fwrite(fid, size(joint_mean{i}{j}, 2), 'int16');
        fwrite(fid, joint_mean{i}{j}, 'double');
        
        fwrite(fid, size(joint_inv_cov{i}{j}, 1), 'int16');
        fwrite(fid, size(joint_inv_cov{i}{j}, 2), 'int16');
        fwrite(fid, joint_inv_cov{i}{j}, 'double');
        
        fclose(fid);
        
        %%
        
        fid = fopen(['joint_normal_', num2str(i-1), '_', num2str(j-1), '.dat'], 'r');
        
        rows = fread(fid, 1, 'int16');
        cols = fread(fid, 1, 'int16');
        joint_mean_ij = fread(fid, rows * cols, 'double');
        joint_mean_ij = reshape(joint_mean_ij, [rows, cols]);
        assert(sum(sum(abs(joint_mean{i}{j} - joint_mean_ij))) == 0);
        
        rows = fread(fid, 1, 'int16');
        cols = fread(fid, 1, 'int16');
        joint_inv_cov_ij = fread(fid, rows * cols, 'double');
        joint_inv_cov_ij = reshape(joint_inv_cov_ij, [rows, cols]);
        assert(sum(sum(abs(joint_inv_cov{i}{j} - joint_inv_cov_ij))) == 0);
        
        fclose(fid);
    end
end
