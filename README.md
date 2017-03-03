# TopologicalDensityEstimation

For details, see https://arxiv.org/abs/1701.09025. Most of the code is to support the evaluation and plots therein. Persons uninterested in reproducing those results can safely focus on tde1d.m and unidec.m. These respectively perform topological density estimation and produce a unimodal decomposition of a function using the algorithm of Baryshnikov and Ghrist; the former calls the latter.

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
