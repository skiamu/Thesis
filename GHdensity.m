function f = GHdensity(param,u,z)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% 1) extract parameters w(k+1) distribution (multivariate)
lambda = param.lambda; Chi = param.Chi; Psi = param.Psi;
mu = param.mu; sigma = param.sigma; gamma = param.gamma;

% 2) compute parameters f(x,u,w(k+1)) distribution (univariate)
mu_bar = u' * mu;
sigma_bar = u' * sigma * u;
gamma_bar =  u' * gamma;

% 3) write explicitly the density function
c = sqrt(Chi*Psi)^(-lambda) * Psi^(lambda) * (Psi+gamma_bar^2 / sigma_bar)^(0.5-lambda) / ...
	(sqrt(2*pi*sigma_bar) * besselk(lambda,sqrt(Chi*Psi)));

factor = sqrt((Chi + (z - mu_bar).^2 / sigma_bar) * (Psi + gamma_bar^2 / sigma_bar));

f = c * besselk(lambda-0.5,factor) .* exp((z - mu_bar)*gamma_bar/sigma_bar) ./ ...
	factor.^(0.5 - lambda);

end

