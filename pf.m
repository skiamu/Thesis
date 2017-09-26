function [ f ] = pf(z,x,u,param)
%pf defines the scalar density function of the rundom variable x(k)*(1+u'*w(k+1))
%   INPUT:
%      z = real variable
%      x = realization portfolio value [scalar]
%      u = vector of asset allocation at time k
%      param = vector of parameters
%   OUTPUT:
%      f = density function
% REMARKS : implement the density also for other models
mu = param.mu; % mean vector
S = param.S; % covariance matrix
f = exp(-0.5*((z - x.*(1 + u'*mu)).^2) ./ (x.^2 * u'*S*u)) ./ ...
	sqrt(2*pi*x.^2*u'*S*u);
end

