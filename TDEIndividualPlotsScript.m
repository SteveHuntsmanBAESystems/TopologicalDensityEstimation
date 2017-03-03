% Script for plotting the densities $f_j$ ($1 \le j \le 6$) and $f_{km}$
% ($1 \le k \le 3$, $1 \le m \le 10$), generating individual sample sets
% from each, performing CV and TDE, and visualizing the results.
%
% NB. The statistics toolbox is required for this script.

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

rng('default');

%% Script configuration
% Linear space
x_0 = -1;
x_1 = 2;
n_x = 500;  % number of evaluation points for PDFs
x = linspace(x_0,x_1,n_x);
% Number of samples
n = 200;
% Kernel
kernelflag = 1;
if kernelflag > 0
    kernelargstr = 'kernel';
    titlestr = 'Gaussian';
elseif kernelflag < 0
    kernelargstr = 'epanechnikov';
    titlestr = 'Epanechnikov';
else
    error('bad kernel flag');
end

%% Generate samples and PDFs
[X,f] = tdepdfsuite(n_x,n);

%% Plot (used to generate f_j.png)
figure;
ha = tight_subplot(5,6,[.01 .01],[.1 .1],[.1 .1]);
for j = 1:6
    axes(ha(j));    
    plot(x,f(j,:),'k');
    xlim([-.5,1.5])
    ymax = 1.1*max(f(j,:));
    ylim([0,ymax]);
    set(ha(j),'XTick',[],'YTick',[]);
end
for j = 7:30, axes(ha(j)); axis off; end

%% Plot (used to generate f_km.png)
figure;
Mvert = 5;
Mhorz = 6;
ha = tight_subplot(Mvert,Mhorz,[.01 .01],[.1 .1],[.1 .1]);
for k = 1:3
    for m = 1:Mhorz
        j = (k-1)*Mhorz+m;
        km = 6+(k-1)*10+m;
        axes(ha(j));    
        plot(x,f(km,:),'k');
        xlim([-.5,1.5]);
        ymax = 1.1*max(f(km,:));
        ylim([0,ymax]);
        set(ha(j),'XTick',[],'YTick',[]);
    end
end

%% Density estimates and plots 
for j = 1:size(X,1)
    %% Optimal kernel density estimate w/ knowledge of PDF
    opt = optimalkernelbandwidth(x,f(j,:),n,kernelflag);
    pd_opt = fitdist(X(j,:)',kernelargstr,'BandWidth',opt.h);
    kde_opt = pdf(pd_opt,x);
    %% Numerically ISE minimizing bandwidth
    temp = vbise(x,f(j,:),X(j,:),kernelflag);
    ind = find(temp.ISE==min(temp.ISE),1,'first');
    h_opt = temp.h(ind);
    pd_numopt = fitdist(X(j,:)',kernelargstr,'BandWidth',h_opt);
    kde_numopt = pdf(pd_numopt,x);
    %% Cross-validation KDE
    cv = cv1d(X(j,:)',kernelflag);
    pd = fitdist(X(j,:)',kernelargstr,'BandWidth',cv.h);
    % THE LINE BELOW IS FOR HISTORICAL REFERENCE ONLY AND DOES NOT DO CV
    % pd = fitdist(X(j,:)',kernelargstr);
    kde = pdf(pd,x);
    %% Topological KDE
    tde = tde1d(X(j,:),kernelflag);
    f_tde = tde.y/(sum(tde.y)*mean(diff(tde.x)));
    %% Plot
    figure;
    plot(x,f(j,:),'k',...
        x,kde_opt,'b',...
        x,kde_numopt,'c',...
        x,kde,'g',...
        tde.x,f_tde,'r');
    xlim([-.5,1.5])
    set(gca,'XTick',[],'YTick',[]);
    title([titlestr,' kernel density estimates of PDF'],...
        'Interpreter','latex');
    legend({'PDF',...
        ['$\hat h_0 =$ ',num2str(opt.h,'%0.4f')],...
        ['$h_{opt} =$ ',num2str(h_opt,'%0.4f')],...
        ['$\hat h_{CV} =$ ',num2str(pd.BandWidth,'%0.4f')],...
        ['$\hat h_{top} =$ ',num2str(tde.h,'%0.4f')]},...
        'Interpreter','latex');
end
