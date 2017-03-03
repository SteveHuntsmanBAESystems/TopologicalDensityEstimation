% Script for plotting performance data for evaluating CV and TDE on the
% densities $f_j$ ($1 \le j \le 6$).

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

figure;

%% Normalization macro
normmacro = ['A = A./repmat(sum(A,1),[size(A,1),1]);',...
    'B = B./repmat(sum(B,1),[size(B,1),1]);'];

%% Plot macro
% "box on" gives inadequate results but retained for tickmarks; manually
% set axis and put in a box since "box on" doesn't give adequate results
plotmacro = ['layertrans2(0:5,[L,Inf],A,B);box on;',... 
    'set(gca,''XTick'',1:5,''XTickLabel'',',...
    '{''25'',''50'',''100'',''200'',''500''});',...
    'hold on;axis([0.5,5.5,lo,hi]);ax = axis;',...
    'line([ax(1),ax(1)]+1e-6,[ax(3),ax(4)],''Color'',''k'');',...
    'line([ax(2),ax(2)]-1e-6,[ax(3),ax(4)],''Color'',''k'');',...
    'line([ax(1),ax(2)],[ax(3),ax(3)]+1e-6,''Color'',''k'');',...
    'line([ax(1),ax(2)],[ax(4),ax(4)]-1e-6,''Color'',''k'');'];

%% c1
for j = 1:6
    % Normalized histogram
    lo = min([min(min(pre_c1_top(:,j,:))),min(min(pre_c1_CV(:,j,:)))]);
    hi = max([max(max(pre_c1_top(:,j,:))),max(max(pre_c1_CV(:,j,:)))]);
    L = linspace(lo,hi,50);
    for k = 1:5 % n = [25,50,100,200,500]
        A(:,k) = histc(squeeze(pre_c1_top(:,j,k)),L,1); 
        B(:,k) = histc(squeeze(pre_c1_CV(:,j,k)),L,1); 
    end
    eval(normmacro);
    % Plot
    subplot(3,6,j);  
    eval(plotmacro);
    title(['$f_',num2str(j),'$'],'Interpreter','latex');
end
subplot(3,6,1);
ylabel('$\hat h - h_{opt}$','Interpreter','latex');

%% c2, c3
for j = 1:6
    % Normalized histogram
    lo = 0;
    hi = max([max(max(pre_c23_top(:,j,:))),max(max(pre_c23_CV(:,j,:)))]);
    L = linspace(lo,hi,50);
    for k = 1:5 % n = [25,50,100,200,500]
        A(:,k) = histc(squeeze(pre_c23_top(:,j,k)),L,1); 
        B(:,k) = histc(squeeze(pre_c23_CV(:,j,k)),L,1); 
    end
    eval(normmacro);
    % Plot
    subplot(3,6,j+6);  
    eval(plotmacro);
end
subplot(3,6,7);
ylabel('ISE$(\hat h)$','Interpreter','latex');

%% c4, c5
for j = 1:6
    % Form histograms (note logspace). NB. The only difference between L1
    % and L2 versions is rescaling    
    lo = -4;
    hi = 0.25;
    L = linspace(lo,hi,50);
    % abs turns out to not have any effect, but included for correctness
    for k = 1:5 % n = [25,50,100,200,500]
        A(:,k) = histc(squeeze(log10(abs(pre_c45_top(:,j,k)))),L,1); 
        B(:,k) = histc(squeeze(log10(abs(pre_c45_CV(:,j,k)))),L,1); 
    end
    eval(normmacro);    
    % Plot
    subplot(3,6,j+12);  
    eval(plotmacro);
    xlabel('$n$','Interpreter','latex');
end
subplot(3,6,13);
ylabel('$c_{45}$','Interpreter','latex');
