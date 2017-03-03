function ISE = ise(x,p,X,h,flag)

% For argument x of a PDF p and samples X drawn from p, this computes the
% integrated square error (ISE) of a kernel density estimate with bandwidth
% h. If flag > 0, the kernel is Gaussian; otherwise, the kernel is
% Epanechnikov.

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

X = X(:);
n = numel(X);

%% Form KDE with bandwidth h
if flag > 0	% Gaussian KDE
    p_hat = zeros(1,numel(x));
    for k = 1:n
        p_hat = p_hat+(1/n)*exp(-.5*((x-X(k))/h).^2)/(h*sqrt(2*pi));
    end
else	% Epanechnikov KDE
    p_hat = zeros(1,numel(x));
    for k = 1:n
        p_hat = p_hat+(1/n)*(.75/h)*max(0,1-((x-X(k))/h).^2);
    end        
end

%% Output
ISE = sum(((p(2:end)-p_hat(2:end)).^2).*diff(x));
