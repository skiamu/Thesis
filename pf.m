function [ f ] = pf(z,x,u,param,model)
%pf defines the scalar density function of the random variable x(k)*(1+u'*w(k+1))
%   INPUT:
%      z = real variable (integration variable) [vector]
%      x = realization portfolio value [scalar]
%      u = vector of asset allocation at time k
%      param = vector of parameters [struct or cell array]
%   OUTPUT:
%      f = density function
% REMARKS:
%      1) when the mixture is used, param is a cell array of the following
%         form {{mu1, Sigma1, lambda1}, ..., {mu_n, Sigma_n, lambda_n}}
%      2) implement the density also for other models

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
	case 'GH' % general case
		% 1) extract parameters w(k+1) distribution (multivariate)
		lambda = param.lambda; Chi = param.Chi; Psi = param.Psi;
		mu = param.mu; sigma = param.sigma; gamma = param.gamma;
		
		% 2) compute parameters f(x,u,w(k+1)) distribution (univariate)
		mu_bar = x * (1 + u' * mu);
		sigma_bar = x^2 * u' * sigma * u;
		gamma_bar = x * u' * gamma;
		
		% 3) write explicitly the density function
		c = sqrt(Chi*Psi)^(-lambda) * Psi^(lambda) * (Psi+gamma_bar^2 / sigma_bar)^(0.5-lambda) / ...
			(sqrt(2*pi*sigma_bar) * besselk(lambda,sqrt(Chi*Psi)));
		
		factor = sqrt((Chi + (z - mu_bar).^2 / sigma_bar) * (Psi + gamma_bar^2 / sigma_bar));
		
		f = c * besselk(lambda-0.5,factor) .* exp((z - mu_bar)*gamma_bar/sigma_bar) ./ ...
			factor.^(0.5 - lambda);
			
	otherwise
		error('invalid model %s',model);
end
end % pf

