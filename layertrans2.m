function h = layertrans2(x,y,A,B)

% Transparency plot of two matrices. This particular code avoids annoyances
% that arise when trying to use surf with transparency options.

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

%% Basic assertions on input sizes
assert(all(size(A)==size(B)),'A and B incompatible');
assert(size(A,2)==numel(x)-1,'x mismatch');
assert(size(A,1)==numel(y)-1,'y mismatch');

%% Affine transformations to unit interval
A = affineunit(A);
B = affineunit(B);

%% Plot
hold on;
for j = 1:size(A,1)
    for k = 1:size(A,2)
        pA = patch([x(k),x(k+1),x(k+1),x(k)]+.5*(x(k+1)-x(k)),...
            [y(j),y(j),y(j+1),y(j+1)]+.5*(y(j+1)-y(j)),'b');
        set(pA,'FaceAlpha',A(j,k),'EdgeAlpha',0);
        pB = patch([x(k),x(k+1),x(k+1),x(k)]+.5*(x(k+1)-x(k)),...
            [y(j),y(j),y(j+1),y(j+1)]+.5*(y(j+1)-y(j)),'r');
        set(pB,'FaceAlpha',B(j,k),'EdgeAlpha',0);
    end
end
xlim(mean(x(1:2))+[min(x),max(x)]);
ylim(mean(y(1:2))+[min(y),max(y)]);

end

%% LOCAL FUNCTION

function A = affineunit(A)

minA = min(A(:));
maxA = max(A(:));
constA = (minA==maxA);  % 1 iff A is constant
A = (A-minA)/(maxA-minA+constA);

end
