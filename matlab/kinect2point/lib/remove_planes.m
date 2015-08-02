function[R, t] = remove_planes(XYZ, param)
B = ransacfitplane(XYZ, param, 1);
B = B ./ norm(B(1:3));

% n = size(XYZ, 2);
% outliers = setdiff(1:n,inliers);

center = -B(4) * B(1:3);
z_axis = B(1:3);
y_axis = normc([0; 1; 0] - z_axis(2) * z_axis);
x_axis = normc(cross(y_axis, z_axis));
R = [x_axis'; y_axis'; z_axis'];
t = -R * center;

pnts = R * XYZ + repmat(t, 1, size(XYZ, 2));
n = size(XYZ, 2);

if(sum((pnts(3, :) < 0)) > n/2)
    disp('Flip...')
    %flip
    R_flip = [1,0,0;0,-1,0;0,0,-1];
%     pnts = R_flip * pnts;
    t = R_flip * t;
    R = R_flip * R;
end
