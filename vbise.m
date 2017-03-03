function y = vbise(x,p,X,flag)

% For argument x of a PDF p and samples X drawn from p, this computes the
% integrated square error (ISE) of kernel density estimates over a suitable
% set of bandwidths. If flag > 0, the kernel is Gaussian; otherwise, the
% kernel is Epanechnikov.

%% Select bandwidths
n = numel(X);
DX = max(X)-min(X);
if DX == 0, error('Dx = 0'); end
h = DX./(1:n);

%% Compute ISE as a function of h
ISE = zeros(1,n);
for j = 1:n
    ISE(j) = ise(x,p,X,h(j),flag);
end

%% Output
y.h = h;
y.ISE = ISE;