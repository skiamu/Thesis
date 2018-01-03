function [ f ] = pfDESext1(z,x,u,J_jump,param)
%pfDES computes the density function of the random variable x(k+1)
%(portdolio value at event number k+1)
%   INPUT:
%      z = indipendent variable
%      x = portfolio value last event
%      u = cash weigth
%      J_jump = jump treshold
%      param = struct of model parameters
%   OUTPUT:
%      f = density computed in z
f = zeros(size(z));

csi = x * J_jump * u;

idx1 = z >= x + csi;

idx2 = z >= x - csi;

mu = param.mu; sigma = param.sigma; r = param.r; % extract parameters

mu_tilde = mu - 0.5 * sigma^2;

p = (exp(2 * mu_tilde * J_jump / sigma^2) - 1 ) / ...
	(2 * sinh(2 * mu_tilde * J_jump / sigma^2)); % probability positive jump

if ~isempty(z(idx1))
	f(idx1) = p * eta((z(idx1) - csi) / x,x,mu_tilde,sigma,r,J_jump);
end

if ~isempty(z(idx2))
	f(idx2) = f(idx2) + (1-p) * eta((z(idx2) + csi) / x,x,mu_tilde,sigma,r,J_jump);
end

f = 2 * cosh(mu_tilde * J_jump / sigma^2) * f;

end % pfDES


function [f] = eta(y,x,mu_tilde,sigma,r,J_jump)
%   INPUT:
%      y = (z+-csi) / x
%      x = ptf value
%      mu_tilde = mu - 0.5 * sigma^2;
%      sigma = volatility
%      r = cash yearly return
%      J_jump = jump size

epsilon = 1e-8; % series truncation tolerance

% t = log(y) / r;
% Ntrunc = sqrt(max(1, -2 * J_jump^2 ./ (pi^2 * sigma^2 * t) .* ...
% 	(log((pi^3 * sigma^2 * t * epsilon) / (4 * J_jump^2))...
% 	- J_jump * mu_tilde / sigma^2))); % number of series terms

% MaxN = min(100,max(Ntrunc));
% N = (1:MaxN);
N = (1:100); % truncate the series at the 100-th terms

% y is a column vector, N must be a row vector. In this case the operation
% y.^N gives a matrix [lenght(y) length(N)]. To get the sum of the partial
% sum we need to sum by rows (e.g. sum(,2))
f = sigma^2 * pi / (4 * J_jump^2) * (sum(N .* (-1).^(N+1) .* (1 / (r * x)) .*  ...
	y.^(-(mu_tilde^2 / (2 * sigma^2) + (sigma^2 * N.^2 * pi^2) / (8 * J_jump^2)) / r - 1)...
	.* sin(pi * N / 2),2));

% the numerical truncation of the series may produce negative values,
% espicially for values of z close do x+csi or x-csi.
f(f<0) = 0;
end % eta



