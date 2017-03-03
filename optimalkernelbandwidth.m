function y = optimalkernelbandwidth(x,f,n,kernelflag)

% For a PDF f(x), this gives the optimal bandwidth for a kernel density
% estimate and the approximate corresponding risk. See Wasserman's All of
% Statistics, Theorem 20.14. If kernelflag > 0, a Gaussian kernel is
% assumed; if kernelflag < 0, an Epanechnikov kernel is assumed; otherwise,
% an error is returned.
%
% WARNING: It is ASSUMED that x is regularly spaced and f is suitably nice.

%% Preliminaries
if numel(x) ~= numel(f)
    error('x and f are incompatibly sized');
else
    f = reshape(f,size(x));
end

%% c1, c2
if kernelflag > 0   % Gaussian
    c1 = 1;
    c2 = .5*sqrt(pi);
elseif kernelflag < 0   % Epanechnikov
    c1 = 1/5;
    c2 = 3/5;
else
    error('bad kernel flag');
end

%% c3
dx = mean(diff(x));
d2fdx2 = conv(f,[1,-2,1],'valid')/dx^2;
c3 = sum((d2fdx2.^2))*dx;

%% Output
y.h = (c1^-.4)*(c2^.2)*(c3^-.2)*(n^-.2);
y.R = .25*(1^4)*(y.h^4)*c3+c2*(y.h^-1)*(n^-1);
