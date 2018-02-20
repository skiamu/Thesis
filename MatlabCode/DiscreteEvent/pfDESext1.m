function [ f ] = pfDESext1(z,x,u,J_jump,param)
%pfDES computes the density function of the random variable x(k+1)
%(portdolio value at event number k+1)
%   INPUT:
%      z = indipendent variable [column vector]
%      x = portfolio value last event [scalar]
%      u = cash weigth [scalar]
%      J_jump = jump treshold
%      param = struct of model parameters (mu, sigma, r and p)
%   OUTPUT:
%      f = density function computed in z
f = zeros(size(z));

csi = x * J_jump * u;

idx1 = z > x + csi + eps;

idx2 = z > x - csi + eps;

mu = param.mu; sigma = param.sigma; r = param.r; % extract parameters
p = param.p;

mu_tilde = mu - 0.5 * sigma^2;

if ~isempty(z(idx1))
	f(idx1) = p * Gamma((z(idx1) - csi) / x,mu_tilde,sigma,r,J_jump);
end

if ~isempty(z(idx2))
	f(idx2) = f(idx2) + (1-p) * Gamma((z(idx2) + csi) / x,mu_tilde,sigma,r,J_jump);
end

f = 2 * cosh(mu_tilde * J_jump / sigma^2) / (r * x) * f;

end % pfDES


function [f] = Gamma(y,mu_tilde,sigma,r,J_jump)
%   INPUT:
%      y = (z+-csi) / x (column vector)
%      x = ptf value
%      mu_tilde = mu - 0.5 * sigma^2;
%      sigma = volatility
%      r = cash yearly return
%      J_jump = jump size
% min_t = 1e-4;
% flag = false;
% epsilon = 1e-8; % series truncation tolerance
% t = min(log(y)/r);
% if t <= min_t
% 	n = length(y); % save the length of y
% 	idx = log(y)/r > min_t; % index of element  not too small
% 	y = y(idx); % remove elemets too small
% 	t = min(log(y)/r); % compute new t
% 	flag = true; % flag there're small elements
% end
% Ntrunc = sqrt(max(1, -8 * J_jump^2 ./ (pi^2 * sigma^2 * t) .* ...
% 	(log((pi^3 * sigma^2 * t * epsilon) / (16 * J_jump^2))...
% 	- J_jump * mu_tilde / sigma^2))); % number of series terms

epsilon = 1e-8;
arg = epsilon * pi * min(y) * log(min(y)) / r;
logArg = log(arg)/log(min(y));
NN = sqrt(-r*8*J_jump^2/(sigma^2*pi^2)*(mu_tilde^2/(2*r*sigma^2) + logArg));
N = (1:NN);
% N = (1:90); % truncate the series at the 100-th terms
global freeMemory
if 8 * NN * length(y) * 1e-6 >= freeMemory % check if there's enogh memory
	f = [0; Gamma(y(2:end),mu_tilde,sigma,r,J_jump)];
else
% y is a column vector, N must be a row vector. In this case the operation
% y.^N gives a matrix [lenght(y) length(N)]. To get the sum of the partial
% sum we need to sum by rows (e.g. sum(,2))
	f = sigma^2 * pi / (4 * J_jump^2) * (sum(N .* (-1).^(N+1) .*  ...
		y.^(-(mu_tilde^2 / (2 * sigma^2) + (sigma^2 * N.^2 * pi^2) / (8 * J_jump^2)) / r - 1)...
		.* sin(pi * N / 2),2));
end
% if flag
% 	g = zeros([n 1]);
% 	g(idx) = f;
% 	f = g;
% end
end % Gamma

