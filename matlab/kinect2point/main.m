close all;
clear all;
addpath 'lib'
addpath 'lib/kovesi'

filename = 'chair003.out';


%% parameters
x_size = 512;
y_size = 424;
n = x_size * y_size;
plane_dist_tol = 50;


%% load data
fileID = fopen(filename);
A = fread(fileID,[x_size y_size],'ushort');
A = A';

max_A = max(A(:));
Depth = double(A) / max_A;
Depth(Depth < 0) = 0;

Depth_filtered = medfilt2(Depth,[8 8]);
% filter_w     = 5;       % bilateral filter half-width
% filter_sigma = [3 0.1]; % bilateral filter standard deviations
% Depth_filtered = bfilter2(Depth_filtered,filter_w,filter_sigma);
A = Depth_filtered * max_A;

depth_image_brightness = 3;
Depth = Depth * depth_image_brightness;
Depth_filtered = Depth_filtered * depth_image_brightness;
hold on; imshow(Depth); title('Depth image'); hold off;
hold on; imshow(Depth_filtered); title('Depth image (filtered)'); hold off;
imwrite(Depth, 'depth_image.png');
imwrite(Depth_filtered, 'depth_image_filtered.png');
fclose(fileID);

%%
XYZ = depth2point(A);

inlier_idx = (XYZ(3, :)' > 0) .* (abs(XYZ(1, :)') < 2000) .* (abs(XYZ(2, :)') < 2000) .* (abs(XYZ(3, :)') < 2000);
inlier_XYZ = XYZ(:, (inlier_idx == 1));
n = sum(inlier_idx == 1);


%% detect ground plane
sample_idx = randi(n, 50000, 1);
sparse_XYZ = inlier_XYZ(:, sample_idx);

[R, t] = remove_planes(sparse_XYZ, plane_dist_tol);
XYZ = R * XYZ + repmat(t, 1, size(XYZ, 2));
inlier_idx(abs(XYZ(3, :)) < plane_dist_tol) = 0;

% s = (1 ./ (max(XYZ(3, :))));
% s = s * 4;
s = 8.4e-04;
XYZ = s * eye(3) * XYZ;

% R_new = vrrotvec2mat([0 0 1 (135 / 180 * pi)]);
% XYZ = R_new * XYZ;
% t = R_new * t;
% R = R_new * R;


%% camera pose
t = s * (-R' * t);
R = R';

R_given = [1 0 0; 0 -1 0; 0 0 -1];
t = R_given * t;
R = R_given * R;

write_camera_pose('camera_pose.txt', R, t);


%% save mesh
px = reshape(XYZ(1, :), y_size, x_size);
py = reshape(XYZ(2, :), y_size, x_size);
pz = reshape(XYZ(3, :), y_size, x_size);
pb = reshape(inlier_idx, y_size, x_size);
v_idx = zeros(y_size, x_size);
count_vertices = 0;

max_edge_length = 0.05;

f = fopen('output.obj', 'w');
for y = 1:y_size
    for x = 1:x_size
        if pb(y, x) == 1
            fprintf(f, 'v %f %f %f\n',...
                px(y, x),...
                py(y, x),...
                pz(y, x));
            
            count_vertices = count_vertices + 1;
            v_idx(y, x) = count_vertices;
        end
    end
end

[grid_x, grid_y] = meshgrid(1:x_size-1, 1:y_size-1);
grid_x = reshape(grid_x, size(grid_x, 1) * size(grid_x, 2), 1);
grid_y = reshape(grid_y, size(grid_y, 1) * size(grid_y, 2), 1);

% v1x v1y v2x v2y v3x v3y
faces_1 = [grid_y, grid_x, grid_y+1, grid_x, grid_y, grid_x+1];
faces_2 = [grid_y+1, grid_x+1, grid_y, grid_x+1, grid_y+1, grid_x];
faces = [faces_1; faces_2];

for i = 1:size(faces, 1)
    face = reshape(faces(i, :), 2, 3)';
    face_vertices = zeros(3, 3);
    
    is_face_valid = true;
    for k = 1:3
        face_vertices(k, :) =...
            [px(face(k, 1), face(k, 2)),...
            py(face(k, 1), face(k, 2)),...
            pz(face(k, 1), face(k, 2))];
    end
    
    for k = 1:3
        k_next = mod(k, 3) + 1;
        if pb(face(k, 1), face(k, 2)) == 0
           is_face_valid = false;
        end
        edge_length = norm(face_vertices(k, :) - face_vertices(k_next, :));
        if edge_length > max_edge_length
            is_face_valid = false;
        end   
    end
    
    if is_face_valid
        fprintf(f, 'f %d %d %d\n',...
            v_idx(face(1, 1), face(1, 2)),...
            v_idx(face(2, 1), face(2, 2)),...
            v_idx(face(3, 1), face(3, 2)));
    end
end

% for y = 1:y_size-1
%     for x = 1:x_size-1
%         if pb(y, x) == 1 && pb(y+1, x) == 1 && pb(y, x+1) == 1
%             fprintf(f, 'f %d %d %d\n',...
%                 v_idx(y, x),...
%                 v_idx(y+1, x),...
%                 v_idx(y, x+1));
%         end
%         
%         if pb(y+1, x+1) == 1 && pb(y, x+1) == 1 && pb(y+1, x) == 1 
%             fprintf(f, 'f %d %d %d\n',...
%                 v_idx(y+1, x+1),...
%                 v_idx(y, x+1),...
%                 v_idx(y+1, x));
%         end
%     end
% end
fclose(f);