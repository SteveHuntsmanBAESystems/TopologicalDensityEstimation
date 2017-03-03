function cv = cv1d(X,kernelflag)

% Cross-validation bandwidth selector for 1D sample data X. Intended for
% comparison with tde1d. If sign(kernelflag) = 1, a Gaussian kernel is
% used; if sign(kernelflag) = -1, an Epanechnikov kernel is used
% 
% Output fields:
%     h,      bandwidth
%     x,      x-values of CV
%     y,      y-values of CV
%     a,      signficance levels
%     l,      lower bounds
%     u,      upper bounds

%% Preliminaries
X = X(:);
n = numel(X);
DX = max(X)-min(X);
if DX == 0
    warning('Dx = 0');
    cv.h = 0;
    cv.x = mean(X);
    cv.y = 1;
    cv.a = 1;
    cv.l = cv.y;
    cv.u = cv.y;
    return;
end

%% Select bandwidths and minimize risk
nh = min(n,100);
cv.h = NaN;
if kernelflag
    %% Select bandwidths
    h = DX./(1:nh);
    L = linspace(min(X),max(X),nh);
    %% Minimize risk
    risk0 = Inf;
    for i = 1:numel(h)  
        risk = cvrisk(X,h(i),kernelflag);
        if risk < risk0
            risk0 = risk;
            cv.h = h(i);
        end
    end
else
    error('kernelflag = 0: use lowriskhist.m instead');
end

%% Output
if kernelflag > 0	% Gaussian
    cv.x = L;   
    cv.y = zeros(1,numel(L));
    for k = 1:n
        cv.y = cv.y+exp(-.5*((L-X(k))/cv.h).^2)/(cv.h*sqrt(2*pi));
    end
    cv.y = cv.y/n;	% normalize
elseif kernelflag < 0	% Epanechnikov
    cv.x = L;   
    cv.y = zeros(1,numel(L));
    for k = 1:n
        cv.y = cv.y+(.75/cv.h)*max(0,1-((L-X(k))/cv.h).^2);
    end
    cv.y = cv.y/n;    % normalize
end

%% Compute (approximate) confidence bands 
% This would benefit from a running variance computation a la Knuth vol. 2,
% p. 232 or en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
d = ceil(log10(n));
cv.a = logspace(-1,-d,d);    % 1-(confidence levels)
cv.l = cell(1,d);  % lower bounds of approximate confidence bands
cv.u = cell(1,d);  % upper bounds of approximate confidence bands
if kernelflag 
    Y = zeros(n,numel(cv.x));
    % See section 20.3 of Wasserman's book All of Statistics 
    if kernelflag > 0	% Gaussian
        omega = 3*cv.h;
        for i = 1:n
            Y(i,:) = exp(-.5*((cv.x-X(i))/cv.h).^2)/(cv.h*sqrt(2*pi));
        end
    elseif kernelflag < 0	% Epanechnikov
        omega = 2*cv.h;
        for i = 1:n
            Y(i,:) = (.75/cv.h)*max(0,1-((cv.x-X(i))/cv.h).^2);
        end
    end
    m = DX/omega;
    q = erfinv((1-cv.a).^(1/m))/sqrt(2);
    s2 = var(Y,0,1);
    se = sqrt(s2/n);
    for j = 1:d
        cv.l{j} = cv.y-q(j)*se;
        cv.u{j} = cv.y+q(j)*se;
    end
else
    % See Theorem 20.10 of AoS and/or Theorem 6.20 of AoNS
    for j = 1:d
        z = erfinv(1-cv.a(j)/cv.h)/sqrt(2);
        c = sqrt(cv.h/n)*z/2;
        cv.l{j} = max(sqrt(cv.y)-c,0).^2;
        cv.u{j} = (sqrt(cv.y)+c).^2;
    end
end