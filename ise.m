function ISE = ise(x,p,X,h,flag)

% For argument x of a PDF p and samples X drawn from p, this computes the
% integrated square error (ISE) of a kernel density estimate with bandwidth
% h. If flag > 0, the kernel is Gaussian; otherwise, the kernel is
% Epanechnikov.

X = X(:);
n = numel(X);

%% Form KDE with bandwidth h
if flag > 0	% Gaussian KDE
    p_hat = zeros(1,numel(x));
    for k = 1:n
        p_hat = p_hat+(1/n)*exp(-.5*((x-X(k))/h).^2)/(h*sqrt(2*pi));
    end
else	% Epanechnikov KDE
    p_hat = zeros(1,numel(x));
    for k = 1:n
        p_hat = p_hat+(1/n)*(.75/h)*max(0,1-((x-X(k))/h).^2);
    end        
end

%% Output
ISE = sum(((p(2:end)-p_hat(2:end)).^2).*diff(x));
