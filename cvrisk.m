function J_hat = cvrisk(X,h,kernelflag)

% Cross-validation estimate of risk for 1D data X and kernel bandwidth h.
% If sign(kernelflag) = 1, a Gaussian kernel is assumed; if
% sign(kernelflag) = -1, an Epanechnikov kernel is assumed; if kernelflag =
% 0, an error suggesting the use of an alternative function is returned.
%
% NB. The exact estimated risk is returned rather than a common (but quite
% accurate) approximation. The runtime penalty for this is negligible.
%
% NB. There are faster ways to compute the risk, at least for Gaussian
% kernels, using gridding/FFT a la Silverman or the fast Gauss transform.
% However, both of these have runtimes that depend on the precision desired
% and introduce a great deal of complexity into the process.

n = numel(X);
if kernelflag ~=0   % kernel
    X = X(:);
    XXh = bsxfun(@minus,X,X')/h;    % XXh(j,k) = (X(j)-X(k))/h
    % K is the kernel; K0 is its value at 0; K2 is the convolution of K
    % with itself
    if sign(kernelflag) == 1    % Gaussian
        XXh2 = XXh.^2;
        K0 = 1/sqrt(2*pi);
        K = exp(-(XXh2)/2)*K0;
        K2 = exp(-(XXh2)/4)*1/sqrt(4*pi);
    else    % Epanechnikov
        aXXh = abs(XXh);
        K0 = 3/4;
        K = (1-aXXh.^2).*(aXXh<=1)*K0;
        K2 = ((2-aXXh).^3).*(aXXh.^2+6*aXXh+4).*(aXXh<=2)*3/160;
    end
    Kn = K2-(2/(n-1))*(n*K-K0);
    J_hat = sum(Kn(:))/(h*n^2);
else    
    error('kernelflag = 0: use lowriskhist.m instead');
end