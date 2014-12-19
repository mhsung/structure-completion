num_cuboids = length(conditional_mean_A);

for i = 1:num_cuboids
    for j = 1:num_cuboids
        if i == j; continue; end;
        fid = fopen(['conditional_normal_', num2str(i-1), '_', num2str(j-1), '.dat'], 'w');
        
        fwrite(fid, size(conditional_mean_A{i}{j}, 1), 'int16');
        fwrite(fid, size(conditional_mean_A{i}{j}, 2), 'int16');
        fwrite(fid, conditional_mean_A{i}{j}, 'double');
        
        fwrite(fid, size(conditional_mean_b{i}{j}, 1), 'int16');
        fwrite(fid, size(conditional_mean_b{i}{j}, 2), 'int16');
        fwrite(fid, conditional_mean_b{i}{j}, 'double');
        
        fwrite(fid, size(conditional_inv_cov{i}{j}, 1), 'int16');
        fwrite(fid, size(conditional_inv_cov{i}{j}, 2), 'int16');
        fwrite(fid, conditional_inv_cov{i}{j}, 'double');
        
        fclose(fid);
        
        %%
        
        fid = fopen(['conditional_normal_', num2str(i-1), '_', num2str(j-1), '.dat'], 'r');
        
        rows = fread(fid, 1, 'int16');
        cols = fread(fid, 1, 'int16');
        conditional_mean_A_ij = fread(fid, rows * cols, 'double');
        conditional_mean_A_ij = reshape(conditional_mean_A_ij, [rows, cols]);
        assert(sum(sum(abs(conditional_mean_A{i}{j} - conditional_mean_A_ij))) == 0);
        
        rows = fread(fid, 1, 'int16');
        cols = fread(fid, 1, 'int16');
        conditional_mean_b_ij = fread(fid, rows * cols, 'double');
        conditional_mean_b_ij = reshape(conditional_mean_b_ij, [rows, cols]);
        assert(sum(sum(abs(conditional_mean_b{i}{j} - conditional_mean_b_ij))) == 0);
        
        rows = fread(fid, 1, 'int16');
        cols = fread(fid, 1, 'int16');
        conditional_inv_cov_ij = fread(fid, rows * cols, 'double');
        conditional_inv_cov_ij = reshape(conditional_inv_cov_ij, [rows, cols]);
        assert(sum(sum(abs(conditional_inv_cov{i}{j} - conditional_inv_cov_ij))) == 0);
        
        fclose(fid);
    end
end
