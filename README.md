# TopologicalDensityEstimation

MATLAB code for topological density estimation. For details, see https://arxiv.org/abs/1701.09025. 

Most of the code here is to support the evaluation and plots in the paper above. Persons uninterested in reproducing those results can safely focus on tde1d.m and unidec.m. These respectively perform topological density estimation and produce a unimodal decomposition of a function using the algorithm of Baryshnikov and Ghrist; the former calls the latter.

The performance data and plot scripts (named "TDE...Script...") call the various functions below and require the statistics toolbox. However, the essential functions tde1d.m and unidec.m (and most of the other functions below) do not require the statistics toolbox.

    cv1d.m is a cross-validation bandwidth selector and calls cvrisk.m
    cvrisk.m is a cross-validation estimate of risk and is called by cv1d.m
    ise.m computes the integrated square error and is called by vbise.m
    layertrans2.m is a transparency plot of two matrices
    lowriskhist.m is a cross-validation histogram function
    optimalkernelbandwidth.m computes the optimal bandwidth for a kernel density estimate
    tde1d.m computes a topological density estimate and calls unidec.m
    tdepdfsuite.m generates samples and PDFs for the scripts
    unidec.m computes a unimodal decomposition and is called by tde1d.m
    vbise.m computes a variable-bandwidth integrated square error and calls ise.m 

Example command sequence:

    %% Evaluate a uniform mixture of m equispaced Gaussians with same variance
    x = linspace(-1,2,500);
    k = 1;
    m = 7;
    mu = (1:m)'/(m+1);
    sigma = reshape(ones(1,m)*(2^-(k+2))*(m+1)^-2,[1,1,m]);
    mix = ones(m,1)/m;
    obj = gmdistribution(mu,sigma,mix);
    f = pdf(obj,x');

    %% Compute unimodal decomposition and perform a rough plot
    u = unidec(f,0); % suppress automatic area plot
    xx = linspace(min(x),max(x),size(u,2));
    figure;
    subplot(1,size(u,1)+1,1);
    area(xx,u'); 
    axis([-.25,1.25,0,1.1*max(f)]);
    for j = 1:size(u,1)
        % Use PP*u instead of u(j,:) to preserve area plot colors
        PP = diag([zeros(1,j-1),1,zeros(1,size(u,1)-j)]);
        subplot(1,size(u,1)+1,j+1);
        area(xx,(PP*u)');
        axis([-.25,1.25,0,1.1*max(f)]);
    end

    %% Perform TDE and plot results
    rng('default');
    X = random(obj,1000)';
    tde = tde1d(X,1); % second arg = 1 for Gaussian kernel
    figure; 
    plot(x,f,'k',tde.x,tde.y,'r');
    legend('PDF','TDE');
    \end{Code}
    We can also superimpose confidence bands:
    \begin{Code}
    hold on;
    for i = 1:numel(tde.a)
        for j = 2:numel(tde.x)
            p = patch(tde.x(j+[-1,0,0,-1]),...
                [tde.l{i}(j+[-1,0]),tde.u{i}(j+[0,-1])],...
                'r','FaceAlpha',2^-(i+1),'EdgeAlpha',0);
        end
    end
    title(['confidence bands: ',regexprep(num2str(1-tde.a),'\s+','; ')]); 

Repeating the above with k = 3 gives a similar but more impressive example.
