function u = unidec(f,plotflag)

% Computes an explicit unimodal decomposition of a (possibly nonnormalized)
% PDF (i.e., a nonnegative function) f with bounded support. For details
% (and some otherwise cryptic notational choices), see Baryshnikov, Y. and
% Ghrist, R. "Unimodal category and topological statistics." NOLTA (2011).

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
if any(f<0), error('f must be nonnegative'); end
f = reshape(f,[1,numel(f)]);
f = [0,f,0];
S = sum(f); 
f = f/S;    % Normalize to unit sum (reverted below)
N = numel(f);
% df = [0,diff(f)];

%% Sweep algorithm
% Step 1
alpha = 1;
u(alpha,:) = zeros(1,N);
g(alpha,:) = f;
% Steps 2/9
while any(g(alpha,:))
    % Step 3: y_alpha = first maximum of g(alpha,:) from left
    y_alpha = find(diff(g(alpha,:))<0,1,'first');
    L = 1:y_alpha;
    R = (y_alpha+1):N;
    % Step 4
    u(alpha+1,L) = g(alpha,L);
    % Step 5 (tweak here courtesy of Jeong-O Jeong)
    df = [0,diff(g(alpha,:))];
    du = min(df(R),0);
    u(alpha+1,R) = u(alpha+1,y_alpha)+cumsum(du);
    % Step 6
    u = max(u,0);
    % Step 7
    alpha = alpha+1;
    % Step 8: alternatively, g(alpha,:) = f-sum(u(1:alpha,:),1)
    g(alpha,:) = g(alpha-1,:)-u(alpha,:);
    g = max(g,0);   % for numerical stability
end

%% Eliminate spurious artifacts
% Revert normalization after elimination
u = S*u(sum(u,2)>eps,:);    
f = S*f;

%% Plot
if plotflag
    figure; 
    area(u');
    hold on
    plot(1:numel(f),f,'k');
    figure; 
    pcolor([[full(log(u));zeros(1,size(u,2))],zeros(size(u,1)+1,1)]);
    shading flat;
    colorbar;
end
