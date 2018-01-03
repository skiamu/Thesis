function [ f ] = pfDESext3(z,x,u,J_jump,param)
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
r0 = 0.05; % initial point
% 2) define mean and std for Integrated O-U r.v.
mu_tilde = @(y) ((r0 - b) * (1 - y.^(a/lambda)) - a * b / lambda * log(y)) / a;
sigma_tilde = @(y) sqrt(sigma^2 / (2 * a^3) * (-2 * a / lambda * log(y) + 4 * y.^(a/lambda) - ...
	y.^(2*a/lambda) - 3));
f1 = zeros(size(z));
f2 = f1;
csi = x * J_jump * u;
idx1 = z >= x + csi;
idx2 = z >= x - csi;

t = (1e-10:0.0001:1-1e-3)';
% 4) compute t-integrals
if ~isempty(z(idx1))
	f1(idx1) = p ./ (z(idx1) - csi) .* trapz(t,Integrand2(t,z(idx1),-1,...
		mu_tilde,sigma_tilde,csi,x));
end
if ~isempty(z(idx2))
    f2(idx2) = (1-p) ./ (z(idx2) + csi) .* trapz(t,Integrand2(t,z(idx2),1,...
		mu_tilde,sigma_tilde,csi,x));
end
f = f1 + f2;
end % pfDES


function I = Integrand2(t,z,sign,mu_tilde,sigma_tilde,csi,x)
assert(iscolumn(t) & isrow(z),'dimension of t or z invalid');
I = normpdf(log((z + sign * csi) / x), mu_tilde(t),...
		sigma_tilde(t)) ;
end
