function Y = multivrandn(u,R,M)
% this function draws M samples from N(u,R)
% where u is the mean vector(row) and R is the covariance matrix which must be positive definite

n = length(u);              % get the dimension
C = chol(R);                % perform cholesky decomp R = C'C
X = randn(M,n);             % draw M samples from N(0,I)

Y = X*C + ones(M,1)*u;      