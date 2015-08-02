function [XYZ] = depth2point(A)

% Kinect v1 parameters.
% fc = [596.23110, 656.18287];
% kc = [-0.01315, 0.00000, -0.00126, 0.00536, 0.00000];
% cc = [375.61455, 437.92128];
% alpha_c = 0.00000;
% height = 480;
% width = 640;

% Kinect v2 parameters.
fc = [377.25386, 373.15983];
cc = [255.25278, 215.64082];
kc = [0.00000, -0.17281, -0.00465, -0.00452, 0.00000];
alpha_c = 0.00000;
height = 424;
width = 512;
err = [0.97735, 1.04026];

dscale = 1;


[u,v] = meshgrid(1:width,1:height);
Z = double(A) * dscale;
%K = [fx 0 cx; 0 fy cy; 0 0 1];
%[xn,yn] = ImageToNormal(u,v,K,rd);

x_kk = [u(:)'; v(:)'];
[x_n] = normalize_pixel(x_kk,fc,cc,kc,alpha_c);
xn = reshape(x_n(1,:)',height,width);
yn = reshape(x_n(2,:)',height,width);

X = xn.* Z;
Y = yn.* Z;

PX = reshape(X, 1, size(X, 1)*size(X, 2));
PY = reshape(Y, 1, size(Y, 1)*size(Y, 2));
PZ = reshape(Z, 1, size(Z, 1)*size(Z, 2));
XYZ = double([PX; PY; PZ]);

% XYZ_norm = sqrt(sum(XYZ .* XYZ, 1));
% idx = XYZ_norm > 1.0E-8;
% XYZ = XYZ(:, idx);
