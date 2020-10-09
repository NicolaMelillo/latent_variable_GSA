function [ mu2, sigma2 ] = getLogNormalValues( what, mu, sigma )

% what = 0
%   mu is the mean OF THE DATA and sigma the standard deviation OF THE DATA
%   (not log transformed). The function will return the parameters of
%   lognormal distribution.
% what = 1
%   mu and sigma are the parameters of the LOGNORMAL DISTRIBUTION (so log
%   transformed). The algorithm will return the
%   mean and the standard deviation of the data (not log transformed).

%%%
% Version 25 October 2018
% by Nicola Melillo
% BMS lab - University of Pavia
% Copyright 2018
%%%

if what==0
    mu2 = log((mu.^2)/sqrt(sigma.^2+mu.^2));
    sigma2 = sqrt(log(sigma.^2./(mu.^2)+1));
elseif what==1
    [mu2,sigma2] = lognstat(mu,sigma);
end

end

