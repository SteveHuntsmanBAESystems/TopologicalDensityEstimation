function y = optimalkernelbandwidth(x,f,n,kernelflag)

% For a PDF f(x), this gives the optimal bandwidth for a kernel density
% estimate and the approximate corresponding risk. See Wasserman's All of
% Statistics, Theorem 20.14. If kernelflag > 0, a Gaussian kernel is
% assumed; if kernelflag < 0, an Epanechnikov kernel is assumed; otherwise,
% an error is returned.
%
% WARNING: It is ASSUMED that x is regularly spaced and f is suitably nice.

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
