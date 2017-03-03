% Script for generating performance data for evaluating CV and TDE on the
% densities $f_j$ ($1 \le j \le 6$) and $f_{km}$ ($1 \le k \le 3$, $1 \le m
% \le 10$).
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

%% Preliminaries and script configuration
% Linear space
x_0 = -1;
x_1 = 2;
n_x = 500;
x = linspace(x_0,x_1,n_x);
% Number of simulation runs
N = 250;
% Number of samples
n = [25,50,100,200,500];	% we omit n = 1000
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

%% Main loop
for aa = 1:numel(n)
    %% Loop over simulation runs
    for ii = 1:N
        disp([aa,ii]);
        %% Generate samples and PDFs
        [X,f] = tdepdfsuite(n_x,n(aa));
        %% Loop over PDFs
        for j = 1:size(X,1) 
            %% Optimal kernel density estimate w/ knowledge of PDF
            opt = optimalkernelbandwidth(x,f(j,:),n(aa),kernelflag);
            if kernelflag > 0
                pd_opt = fitdist(X(j,:)',kernelargstr,'BandWidth',opt.h);
            elseif kernelflag < 0
                pd_opt = fitdist(X(j,:)','Kernel','Kernel',kernelargstr,...
                    'BandWidth',opt.h);
            end
            kde_opt = pdf(pd_opt,x);
            %% Numerically ISE minimizing bandwidth
            temp = vbise(x,f(j,:),X(j,:),kernelflag);
            ind = find(temp.ISE==min(temp.ISE),1,'first');
            h_opt = temp.h(ind);
            if kernelflag > 0
                pd_numopt = fitdist(X(j,:)',kernelargstr,'BandWidth',h_opt);
            elseif kernelflag < 0
                pd_numopt = fitdist(X(j,:)','Kernel','Kernel',...
                    kernelargstr,'BandWidth',h_opt);
            end
            kde_numopt = pdf(pd_numopt,x);
            %% Cross-validation KDE
            cv = cv1d(X(j,:)',kernelflag);
            if kernelflag > 0
                % % For MATLAB estimator
                % pd = fitdist(X(j,:)',kernelargstr);   
                pd = fitdist(X(j,:)',kernelargstr,'BandWidth',cv.h);
            elseif kernelflag < 0
                % % For MATLAB estimator
                % pd = fitdist(X(j,:)','Kernel','Kernel',kernelargstr);
                pd = fitdist(X(j,:)','Kernel','Kernel',kernelargstr,'BandWidth',cv.h);
            end
            kde = pdf(pd,x);
            h_hat = pd.BandWidth;
            pre_c1_CV(ii,j,aa) = h_hat-h_opt;
            pre_c23_CV(ii,j,aa) = ise(x,f(j,:),X(j,:),h_hat,kernelflag);
            pre_c45_CV(ii,j,aa) = ...
                pre_c23_CV(ii,j,aa)-ise(x,f(j,:),X(j,:),h_opt,kernelflag);
            %% Topological KDE
            tde = tde1d(X(j,:),kernelflag);
            f_tde = tde.y/(sum(tde.y)*mean(diff(tde.x)));
            h_hat = tde.h;
            pre_c1_top(ii,j,aa) = h_hat-h_opt;
            pre_c23_top(ii,j,aa) = ise(x,f(j,:),X(j,:),h_hat,kernelflag);
            pre_c45_top(ii,j,aa) = ...
                pre_c23_top(ii,j,aa)-ise(x,f(j,:),X(j,:),h_opt,kernelflag);
            %% Unimodal category and number of local maxima
            ucat_CV(ii,j,aa) = nnz(sum(unidec(kde,0),2)>sqrt(eps));
            lmax_CV(ii,j,aa) = nnz(diff(sign(diff([0,kde,0])))<0);
            ucat_top(ii,j,aa) = tde.mfuc;
            lmax_top(ii,j,aa) = nnz(diff(sign(diff([0,f_tde,0])))<0);
        end  
    end
end
