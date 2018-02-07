function [ f ] = pfDESext2(z,x,u,J_jump,param)
%pfDESext2 computes the density function of the random variable x(k+1)
%(portdolio value at event number k+1)
%   INPUT:
%      z = indipendent variable
%      x = portfolio value last event
%      u = cash weigth
%      J_jump = jump treshold
%      param = struct of model parameters
%   OUTPUT:
%      f = density computed in z

% 1) extract parameters
a = param.a; b = param.b; sigma = param.sigma; % vasicek parameters
p = param.p; lambda = param.lambda; % basic parameters
r0 = 0.055; % initial point
% 2) define mean and std for Integrated O-U r.v.
mu_tilde = @(t) ((r0 - b) * (1 - exp(-a * t)) + a * b * t) / a;
sigma_tilde = @(t) sqrt(sigma^2 / (2 * a^3) * (2 * a * t + 4 * exp(-a * t) - ...
	exp(-2 * a * t) - 3));

% 3) tau density
f_tau = @(t) lambda * exp(-lambda * t);

f = zeros(size(z));
csi = x * J_jump * u;
idx1 = z >= x + csi;
idx2 = z >= x - csi;

t = (-2:0.001/2:2)';
% 4) compute t-integrals
if ~isempty(z(idx1))
	psi = @(t,z) Integrand(exp(t - exp(-t)),z,f_tau,-1,mu_tilde,sigma_tilde,csi,x) .* ...
		exp(t - exp(-t)) .* (1+exp(-t));
	f(idx1) = p ./ (z(idx1) - csi) .* trapz(t,psi(t,z(idx1)));
end
if ~isempty(z(idx2))
	psi = @(t,z) Integrand(exp(t - exp(-t)),z,f_tau,1,mu_tilde,sigma_tilde,csi,x) .* ...
		exp(t - exp(-t)) .* (1+exp(-t));
	f(idx2) = f(idx2) + (1-p) ./ (z(idx2) + csi) .* trapz(t,psi(t,z(idx2)));
end
end % pfDES


function I = Integrand(t,z,f_tau,sign,mu_tilde,sigma_tilde,csi,x)
sigma_tilde = sigma_tilde(t);
t = t(sigma_tilde>0); % keep only times where sigma_tilde is not zero
I =  sparse((exp(-0.5 * ((log((z + sign * csi) / x) - mu_tilde(t)) ./ sigma_tilde).^2 )...
	./ (sqrt(2*pi) * sigma_tilde)  ) .* (f_tau(t)));
end
% * ones([1 length(z)]))
