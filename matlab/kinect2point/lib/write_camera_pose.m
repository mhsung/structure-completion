function [] = write_camera_pose(filename, R, t)
    fid = fopen(filename, 'w');
    for i = 1:3
        for j = 1:3
            fprintf(fid, '%f\n', R(j, i));
        end
        fprintf(fid, '%f\n', 0);
    end
    for i = 1:3
        fprintf(fid, '%f\n', t(i));
    end
    fprintf(fid, '%f\n', 1);
    fclose(fid);