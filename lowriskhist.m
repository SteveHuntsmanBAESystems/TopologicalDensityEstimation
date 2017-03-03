function H = lowriskhist(x)

% Histogram with bins nonparametrically determined by cross-validation
% estimator / estimated risk. See 20.2 of All of Statistics and/or 6.2 of
% All of Nonparametric Statistics, both by Wasserman. Note that (20.14) in
% AoS has a typo: it is missing an h in the denominator of the second term,
% which is included here (and is already in the corresponding (6.16) of
% AoNS).
% 
% Input:	x       A vector of data to be binned
%
% Output:   H       A struct with the following fields:
%           H.x     The optimal bin centers
%           H.y     The corresponding counts
%           H.a     1-(confidence levels)
%           H.l     Lower bounds of approximate confidence bands
%           H.u     Upper bounds of approximate confidence bands
%
% Note that this function can be quite slow, since it has to compute a
% linear number of histograms and therefore requires overall quadratic
% time. For large data sets, running this on subsamples and
% cross-validating the results would seem to be the right thing to do. 

n = numel(x);

%% Compute cross-validation estimator / estimated risk
% While the number of bins should notionally be O(n^(1/3)), the hidden
% constant is not really knowable, so we use m_max = n vs. (e.g.) sqrt(n).
% Note that the example in Wasserman illustrates the need for this
% generality.
m_max = n;
J_hat = zeros(1,m_max);             % cross-validation estimator of risk
for m = 1:m_max
    L = linspace(min(x),max(x),m);	% bins
    h = 1/m;                        % bandwidth
    p_hat = histc(x,L)/n;
    J_hat(m) = 2/((n-1)*h)-((n+1)/((n-1)*h))*sum(p_hat.^2);
end

%% Perform bias-variance tradeoff
m_opt = find(J_hat==min(J_hat),1,'first');
H.x = linspace(min(x),max(x),m_opt);    % optimal bin centers
H.y = histc(x,H.x);                     % corresponding counts

%% Compute confidence bands 
% See Theorem 20.10 of AoS and/or Theorem 6.20 of AoNS
d = ceil(log10(n));
H.a = logspace(-1,-d,d);    % 1-(confidence levels)
H.l = cell(1,d);  % lower bounds of approximate confidence bands
H.u = cell(1,d);  % upper bounds of approximate confidence bands
for j = 1:d
    z = sqrt(2)*erfinv(1-H.a(j)/m_opt);
    c = sqrt(m_opt/n)*z/2;
    H.l{j} = max(sqrt(H.y)-c,0).^2;
    H.u{j} = (sqrt(H.y)+c).^2;
end
