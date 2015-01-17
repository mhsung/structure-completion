function[inv_cov] = pca_inv_cov(cov)
[V,D,E] = pcacov(cov);
pcacov_range = find(cumsum(E) <= 99.9);
V = V(:, pcacov_range);
D = D(pcacov_range);

% Covariance of (X * V) = diag(D)
assert(sum(sum(abs((V'*cov*V - diag(D))))) < 0.001);

% Mahalanobis distance in (X * V):
% (X * V) * inv(diag(D)) * (X * V)'
% X * (V * inv(diag(D)) * V') * X'
inv_cov = V * diag(1.0 ./ D) * V';