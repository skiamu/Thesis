function [ f ] = pf(z,x,u,param,model)
%pf defines the scalar density function of the random variable x(k)*(1+u'*w(k+1))
%   INPUT:
%      z = real variable
%      x = realization portfolio value [scalar]
%      u = vector of asset allocation at time k
%      param = vector of parameters
%   OUTPUT:
%      f = density function
% REMARKS : implement the density also for other models

switch model
	case 'Gaussian'
		mu = param.mu; % mean vector
		S = param.S; % covariance matrix
		f = exp(-0.5*((z - x.*(1 + u'*mu)).^2) ./ (x.^2 * u'*S*u)) ./ ...
			sqrt(2*pi*x.^2*u'*S*u);
	case 'Mixture'
		n = length(param); % param is as long as the mixture addendum
		f = 0;
		for i = 1 : n
			mu = param{i}{1};
			S = param{i}{2};
			lambda = param{i}{3};
			f = f + lambda * exp(-0.5*((z - x.*(1 + u'*mu)).^2) ./ (x.^2 * u'*S*u)) ./ ...
				sqrt(2*pi*x.^2*u'*S*u);
		end
	
	case 'Levy'
		
	otherwise
		error('invalid model %s',model);
end
end % pf

