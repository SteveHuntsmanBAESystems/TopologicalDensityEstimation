function tde = tde1d(X,kernelflag)

% One-dimensional topological density estimate of sample data X. If
% kernelflag = 0, a histogram is returned; if sign(kernelflag) = 1, a
% Gaussian kernel density estimate is returned; if sign(kernelflag) = -1,
% an Epanechnikov kernel density estimate is returned.
% 
% Output fields:
%     h,      bandwidth
%     x,      x-values of TDE
%     y,      y-values of TDE
%     uc,     unimodal category
%     mfuc,   most frequent unimodal category
%     a,      signficance levels
%     l,      lower bounds
%     u,      upper bounds

% NB. The use of histc vs. histcounts is deliberate to accomodate legacy
% installations of MATLAB.

%% Preliminaries
X = X(:);
n = numel(X);
DX = max(X)-min(X);
if DX == 0
    warning('Dx = 0');
    tde.h = 0;
    tde.x = mean(X);
    tde.y = 1;
    tde.uc = 1;
    tde.mfuc = 1;
    tde.a = 1;
    tde.l = tde.y;
    tde.u = tde.y;
    return;
end

%% Select bandwidths/bin numbers
nh = min(n,100);
if kernelflag
    %% Select bandwidths
    h = DX./(1:nh);
else
    %% Select bin numbers
    h = 1:nh;
end
L = linspace(min(X),max(X),nh);

%% Compute unimodal category as a function of h
uc = zeros(1,nh);
for j = 1:nh
    %% Form density estimate corresponding to h(j)
    if kernelflag > 0	% Gaussian KDE
        p_hat = zeros(1,numel(L));
        for k = 1:n
            p_hat = p_hat+exp(-.5*((L-X(k))/h(j)).^2)/(h(j)*sqrt(2*pi));
        end
        p_hat = p_hat/n;	% normalize
    elseif kernelflag < 0 % Epanechnikov KDE
        p_hat = zeros(1,numel(L));
        for k = 1:n
            p_hat = p_hat+(.75/h(j))*max(0,1-((L-X(k))/h(j)).^2);
        end   
        p_hat = p_hat/n;	% normalize
    else	% histogram
        p_hat = histc(X,linspace(min(L),max(L),h(j)))/n;
    end        
    %% Compute unimodal category
    u = unidec(p_hat,0);
    uc(j) = numel(sum(u,2)>0);
end

%% Find most frequent unimodal category
uuc = unique(uc);
fuc = zeros(1,numel(uuc));
for j = 1:numel(uuc)
    fuc(j) = nnz(uc==uuc(j));
end
mfuc = uuc(find(fuc==max(fuc),1,'first'));

%% Find central number of bins for the most frequent unimodal category
temp = find(uc==mfuc);
h_opt = h(temp(ceil(numel(temp)/2)));
tde.h = h_opt;

%% Output
if kernelflag > 0	% Gaussian
    tde.x = L;   
    tde.y = zeros(1,numel(L));
    for k = 1:n
        tde.y = tde.y+exp(-.5*((L-X(k))/h_opt).^2)/(h_opt*sqrt(2*pi));
    end
    tde.y = tde.y/n;	% normalize
elseif kernelflag < 0	% Epanechnikov
    tde.x = L;   
    tde.y = zeros(1,numel(L));
    for k = 1:n
        tde.y = tde.y+(.75/h_opt)*max(0,1-((L-X(k))/h_opt).^2);
    end
    tde.y = tde.y/n;    % normalize
else    % histogram
    tde.x = linspace(min(L),max(L),h_opt);	% optimal bin centers
    tde.y = histc(X,tde.x)/n;               % corresp. normalized counts
end
tde.uc = uc;
tde.mfuc = mfuc;

%% Compute (approximate) confidence bands 
% This would benefit from a running variance computation a la Knuth vol. 2,
% p. 232 or en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
d = ceil(log10(n));
tde.a = logspace(-1,-d,d);    % 1-(confidence levels)
tde.l = cell(1,d);  % lower bounds of approximate confidence bands
tde.u = cell(1,d);  % upper bounds of approximate confidence bands
if kernelflag 
    Y = zeros(n,numel(tde.x));
    % See section 20.3 of Wasserman's book All of Statistics 
    if kernelflag > 0	% Gaussian
        omega = 3*tde.h;
        for i = 1:n
            Y(i,:) = exp(-.5*((tde.x-X(i))/tde.h).^2)/(tde.h*sqrt(2*pi));
        end
    elseif kernelflag < 0	% Epanechnikov
        omega = 2*tde.h;
        for i = 1:n
            Y(i,:) = (.75/tde.h)*max(0,1-((tde.x-X(i))/tde.h).^2);
        end
    end
    m = DX/omega;
    q = erfinv((1-tde.a).^(1/m))/sqrt(2);
    s2 = var(Y,0,1);
    se = sqrt(s2/n);
    for j = 1:d
        tde.l{j} = tde.y-q(j)*se;
        tde.u{j} = tde.y+q(j)*se;
    end
else
    % See Theorem 20.10 of AoS and/or Theorem 6.20 of AoNS
    for j = 1:d
        z = erfinv(1-tde.a(j)/h_opt)/sqrt(2);
        c = sqrt(h_opt/n)*z/2;
        tde.l{j} = max(sqrt(tde.y)-c,0).^2;
        tde.u{j} = (sqrt(tde.y)+c).^2;
    end
end