function [X,f] = tdepdfsuite(n_x,n)

% Generate samples X and PDFs f for evaluation of topological density
% estimation. n_x is the number of points to evaluate PDFs at; n is the
% number of samples for each PDF.
%
% NB. The statistics toolbox is required for this function.
%
% Example:
%   [X,f] = tdepdfsuite(500,200);
 
%% Preliminaries
% Linear space
x_0 = -1;
x_1 = 2;
x = linspace(x_0,x_1,n_x);
% Samples
X = [];

%% Laplace distribution $f_1$
j = size(X,1)+1;
% Parameters
mu = .5;
b = .125;
theta{j} = [mu,b];
% Sample data
obj = makedist('Uniform');
U = random(obj,1,n)-.5;
X(j,:) = theta{j}(1)-theta{j}(2)*sign(U).*log(1-2*abs(U));
% Explicit PDF
f(j,:) = .5*exp(-abs(x-theta{j}(1))/theta{j}(2))/theta{j}(2);

%% Simple gamma distribution $f_2$
j = size(X,1)+1;
% Parameters
b = 1.5;
a = b^2;    
ell = 5;
shape = a;
scale = 1/(ell*b);
theta{j} = [shape,scale];
% Sample data
obj = makedist('Gamma','a',shape,'b',scale);
X(j,:) = random(obj,1,n);
% Explicit PDF
f(j,:) = pdf(obj,x);

%% Mixture of three gamma distributions $f_3$
j = size(X,1)+1;
% Parameters
mix = ones(3,1)/3;
cmix = [0;cumsum(mix(1:end-1))];
b = [1.5;3;6];
a = b.^2;
ell = 8;
shape = a;
scale = 1./(ell*b);
theta{j} = [mix,shape,scale];
% Sample data
for i = 1:3
    obj(i) = makedist('Gamma','a',shape(i),'b',scale(i));
end
X(j,:) = zeros(1,n);
for k = 1:n
    i = find(cmix<rand,1,'last');
    X(j,k) = random(obj(i),1,1);
end 
% Explicit PDF
f(j,:) = zeros(1,n_x);
for i = 1:3
    f(j,:) = f(j,:)+mix(i)*pdf(obj(i),x);
end

%% Simple normal distribution $f_4$
j = size(X,1)+1;
% Parameters
mu = 0.5;
sigma = 0.2;
theta{j} = [mu,sigma];
% Sample data
obj = gmdistribution(mu,sigma'*sigma);
X(j,:) = random(obj,n)';
% Explicit PDF
f(j,:) = pdf(obj,x');

%% Mixture of two normal distributions $f_5$
j = size(X,1)+1;
% Parameters
mu = [0.35;0.65];
sigma = cat(3,0.1^2,0.1^2);
mix = ones(2,1)/2;
theta{j} = [mix,mu,[sigma(1);sigma(2)]];
% Sample data and explicit PDF
obj = gmdistribution(mu,sigma,mix);
X(j,:) = random(obj,n)';
% Explicit PDF
f(j,:) = pdf(obj,x');

%% Mixture of three normal distributions $f_6$
j = size(X,1)+1;
% Parameters
mu = [0.25;0.5;0.75];
sigma = cat(3,0.075^2,0.075^2,0.075^2);
mix = ones(3,1)/3;
theta{j} = [mix,mu,[sigma(1);sigma(2);sigma(3)]];
% Sample data and explicit PDF
obj = gmdistribution(mu,sigma,mix);
X(j,:) = random(obj,n)';
% Explicit PDF
f(j,:) = pdf(obj,x');

%% Mixtures of varying numbers of normal distributions $f_{km}$
for k = 1:3
    for m = 1:10
        j = size(X,1)+1;
        % Parameters
        mu = (1:m)'/(m+1);
        s2 = (2^-(k+2))*(m+1)^-2;
        sigma = reshape(s2*ones(1,m),[1,1,m]);
        mix = ones(m,1)/m;
        theta{j} = [mix,mu,reshape(sigma,[m,1])];
        % Sample data and explicit PDF
        obj = gmdistribution(mu,sigma,mix);
        X(j,:) = random(obj,n)';
        % Explicit PDF
        f(j,:) = pdf(obj,x');       
    end
end
