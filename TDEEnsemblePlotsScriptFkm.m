% Script for plotting performance data for evaluating CV and TDE on the
% densities $f_{km}$ ($1 \le k \le 3$, $1 \le m \le 10$).

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

%% Normalization macro
normmacro = ['A = A./repmat(sum(A,1),[size(A,1),1]);',...
    'B = B./repmat(sum(B,1),[size(B,1),1]);'];

%% Plot macro
% "box on" gives inadequate results but retained for tickmarks; manually
% set axis and put in a box since "box on" doesn't give adequate results
plotmacro = ['layertrans2(0:10,[L,Inf],A(:,ind),B(:,ind));box on;',... 
    'set(gca,''XTick'',1:10,''XTickLabel'',',...
    '{''1'','''','''','''','''','''','''','''','''',''10''});',...
    'hold on;axis([0.5,10.5,lo,hi]);ax = axis;',...
    'line([ax(1),ax(1)]+1e-6,[ax(3),ax(4)],''Color'',''k'');',...
    'line([ax(2),ax(2)]-1e-6,[ax(3),ax(4)],''Color'',''k'');',...
    'line([ax(1),ax(2)],[ax(3),ax(3)]+1e-6,''Color'',''k'');',...
    'line([ax(1),ax(2)],[ax(4),ax(4)]-1e-6,''Color'',''k'');'];

%% Loop over values of n
for j = 1:5 % n = [25,50,100,200,500]
    figure;

    %% c1
    % Form histograms
    lo = min([min(min(pre_c1_top(:,:,j))),min(min(pre_c1_CV(:,:,j)))]);
    hi = max([max(max(pre_c1_top(:,:,j))),max(max(pre_c1_CV(:,:,j)))]);
    L = linspace(lo,hi,50);
    A = histc(pre_c1_top(:,:,j),L,1); 
    B = histc(pre_c1_CV(:,:,j),L,1); 
    eval(normmacro);
    % Transparency plots
    subplot(3,3,1); ind = 7:16; eval(plotmacro);
    ylabel('$k=1$','Interpreter','latex');
    title('$\hat h - h_{opt}$','Interpreter','latex');
    subplot(3,3,4); ind = 17:26; eval(plotmacro);
    ylabel('$k=2$','Interpreter','latex');
    subplot(3,3,7); ind = 27:36; eval(plotmacro);
    xlabel('$m$','Interpreter','latex');
    ylabel('$k=3$','Interpreter','latex');

    %% c2, c3
    % Form histograms
    lo = 0;
    hi = max([max(max(pre_c23_top(:,:,j))),max(max(pre_c23_CV(:,:,j)))]);
    L = linspace(lo,hi,50);
    A = histc(pre_c23_top(:,:,j),L,1); 
    B = histc(pre_c23_CV(:,:,j),L,1); 
    eval(normmacro);
    % Transparency plots
    subplot(3,3,2); ind = 7:16; eval(plotmacro);
    title('ISE$(\hat h)$','Interpreter','latex');
    subplot(3,3,5); ind = 17:26; eval(plotmacro);
    subplot(3,3,8); ind = 27:36; eval(plotmacro);
    xlabel('$m$','Interpreter','latex');

    %% c4, c5
    % Form histograms (note logspace). NB. The only difference between L1
    % and L2 versions is rescaling
    % The following is roughly equivalent in practice to 
    %   temp = cat(3,pre_c45_top(:,:,j),pre_c45_CV(:,:,j));
    %   qlog = .05; % quantile margin for log cutoff
    %   lo = log10(quantile(temp(:),qlog));
    %   hi = max([0,log10(max(temp(:)))]);
    % but has the distinct advantage of being constant
    lo = -4;    
    hi = 0.25;
    L = linspace(lo,hi,50);
    % abs turns out to not have any effect, but included for correctness
    A = histc(log10(abs(pre_c45_top(:,:,j))),L,1); 
    B = histc(log10(abs(pre_c45_CV(:,:,j))),L,1); 
    eval(normmacro);
    % Transparency plots
    subplot(3,3,3); ind = 7:16; eval(plotmacro);
	title('$c_{45}$','Interpreter','latex');
    subplot(3,3,6); ind = 17:26; eval(plotmacro);
    subplot(3,3,9); ind = 27:36; eval(plotmacro);
    xlabel('$m$','Interpreter','latex');
    
    %% Figure output
    print('-dpng',[titlestr,'_',num2str(j),'.png']);
end

%% ucat
for i = 1:size(f,1)
    u = unidec(f(i,:),0); 
    ucat(i) = size(u,1); 
end
for j = 1:5
    figure;
    % Form histograms
    lo = 1;
    hi = 21;
    L = 0.5:1:20.5;
    A = histc(ucat_top(:,:,j),L,1); 
    B = zeros(size(A)); B(end,:) = 1; 
    eval(normmacro);    
    % Transparency plots
    subplot(3,3,1); ind = 7:16; eval(plotmacro);
    ylabel('$k=1$','Interpreter','latex');
	title('ucat$_{TDE}$','Interpreter','latex');
    subplot(3,3,4); ind = 17:26; eval(plotmacro);
    ylabel('$k=2$','Interpreter','latex');
    subplot(3,3,7); ind = 27:36; eval(plotmacro);    
    xlabel('$m$','Interpreter','latex');
    ylabel('$k=3$','Interpreter','latex');
    A = B;
    B = histc(ucat_CV(:,:,j),L,1); 
    eval(normmacro);    
    subplot(3,3,2); ind = 7:16; eval(plotmacro);
	title('ucat$_{CV}$','Interpreter','latex');
    subplot(3,3,5); ind = 17:26; eval(plotmacro);
    subplot(3,3,8); ind = 27:36; eval(plotmacro);    
    xlabel('$m$','Interpreter','latex');
    % Diagonal plots
    A = histc(ucat_top(:,:,j),L,1); 
    B = histc(ucat_CV(:,:,j),L,1);
    subplot(3,3,3); ind = 7:16; 
    for i = 1:numel(ind)
        tempA(i) = A(ucat(ind(i)),ind(i))/N;
        tempB(i) = B(ucat(ind(i)),ind(i))/N;
    end
    plot(1:10,tempA,'bo',1:10,tempB,'rx'); xlim([0,11]);
    set(gca,'XTick',1:10,'XTickLabel',{'1','','','','','','','','','10'});
	title('Pr$($ucat$(\hat f)=$ucat$(f))$','Interpreter','latex');
    subplot(3,3,6); ind = 17:26; 
    for i = 1:numel(ind)
        tempA(i) = A(ucat(ind(i)),ind(i))/N;
        tempB(i) = B(ucat(ind(i)),ind(i))/N;
    end
    plot(1:10,tempA,'bo',1:10,tempB,'rx'); xlim([0,11]);
    set(gca,'XTick',1:10,'XTickLabel',{'1','','','','','','','','','10'});
    subplot(3,3,9); ind = 27:36; 
    for i = 1:numel(ind)
        tempA(i) = A(ucat(ind(i)),ind(i))/N;
        tempB(i) = B(ucat(ind(i)),ind(i))/N;
    end
    plot(1:10,tempA,'bo',1:10,tempB,'rx'); xlim([0,11]);
    set(gca,'XTick',1:10,'XTickLabel',{'1','','','','','','','','','10'});
    xlabel('$m$','Interpreter','latex');
    %% Figure output
    print('-dpng',['UCAT',titlestr,'_',num2str(j),'.png']);
end

%% lmax
for j = 1:5
    figure;
    % Form histograms
    lo = 1;
    hi = 21;
    L = 0.5:1:20.5;
    A = histc(lmax_top(:,:,j),L,1); 
    B = zeros(size(A)); B(end,:) = 1; 
    eval(normmacro);    
    % Transparency plots
    subplot(3,3,1); ind = 7:16; eval(plotmacro);
    ylabel('$k=1$','Interpreter','latex');
	title('$\#_{max,TDE}$','Interpreter','latex');
    subplot(3,3,4); ind = 17:26; eval(plotmacro);
    ylabel('$k=2$','Interpreter','latex');
    subplot(3,3,7); ind = 27:36; eval(plotmacro);    
    xlabel('$m$','Interpreter','latex');
    ylabel('$k=3$','Interpreter','latex');
    A = B;
    B = histc(lmax_CV(:,:,j),L,1); 
    eval(normmacro);    
    subplot(3,3,2); ind = 7:16; eval(plotmacro);
	title('$\#_{max,CV}$','Interpreter','latex');
    subplot(3,3,5); ind = 17:26; eval(plotmacro);
    subplot(3,3,8); ind = 27:36; eval(plotmacro);    
    xlabel('$m$','Interpreter','latex');
    % Diagonal plots
    A = histc(lmax_top(:,:,j),L,1); 
    B = histc(lmax_CV(:,:,j),L,1);
    subplot(3,3,3); ind = 7:16; 
    tempA = diag(A(:,ind))/N; tempB = diag(B(:,ind))/N;
    plot(1:10,tempA,'bo',1:10,tempB,'rx'); xlim([0,11]);
    set(gca,'XTick',1:10,'XTickLabel',{'1','','','','','','','','','10'});
	title('Pr$(\#_{max} = m)$','Interpreter','latex');
    subplot(3,3,6); ind = 17:26; 
    tempA = diag(A(:,ind))/N; tempB = diag(B(:,ind))/N;
    plot(1:10,tempA,'bo',1:10,tempB,'rx'); xlim([0,11]);
    set(gca,'XTick',1:10,'XTickLabel',{'1','','','','','','','','','10'});
    subplot(3,3,9); ind = 27:36; 
    tempA = diag(A(:,ind))/N; tempB = diag(B(:,ind))/N;
    plot(1:10,tempA,'bo',1:10,tempB,'rx'); xlim([0,11]);
    set(gca,'XTick',1:10,'XTickLabel',{'1','','','','','','','','','10'});
    xlabel('$m$','Interpreter','latex');
    %% Figure output
    print('-dpng',['LMAX',titlestr,'_',num2str(j),'.png']);
end
