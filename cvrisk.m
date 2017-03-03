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

% Copyright (c) 2017, BAE Systems. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright notice,
% this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright
% notice,
%   this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
% 
% * Neither the name of the copyright holder nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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