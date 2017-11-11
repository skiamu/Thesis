function f = pf(z,x,u,param,model)
%pf defines the scalar density function of the random variable x(k)*(1+u'*w(k+1))
%   INPUT:
%      z = real variable (integration variable) [vector]
%      x = realization portfolio value [scalar]
%      u = vector of asset allocation at time k [vector]
%      param = vector of parameters [struct array]
%   OUTPUT:
%      f = density function in z

switch model
	case 'Gaussian'
		mu = x * (1 + u' * param.mu);
		sigma = sqrt(x^2 * u' * param.S * u);
		f = exp(-0.5 * (z - mu).^2 / sigma^2) / sqrt(2 * pi * sigma^2);
	case 'Mixture'
		n = length(param);
		f = 0;
		for i = 1 : n
			mu = x * (1 + u' * param(i).mu);
			sigma = sqrt(x^2 * (u' * param(i).S * u));
			f = f + param(i).lambda * normpdf(z,mu,sigma);
		end
	case {'GH','H','NIG','t'}
		
		if strcmp(model,'t')
			% 1) extract parameters w(k+1) distribution (multivariate)
			nu = param.nu; 
			mu = param.mu; sigma = param.sigma; gamma = param.gamma;
			
			% 2) compute parameters f(x,u,w(k+1)) distribution (univariate)
			mu_bar = x * (1 + u' * mu);
			sigma_bar = x^2 * u' * sigma * u;
			gamma_bar = x * u' * gamma;
			
			% 3) write explicitly the density function
			c = ((nu-2)^(nu/2) * (gamma_bar^2/sigma_bar)^((nu+1)/2)) /...
				(sqrt(2*pi*sigma_bar) * feval('gamma',(nu/2)) * 2^(nu/2 - 1));
			
			factor = sqrt((nu-2 + (z - mu_bar).^2 / sigma_bar) * (gamma_bar^2/sigma_bar));
			
			f = c * besselk((nu+1)/2,factor) .* exp((z - mu_bar)*gamma_bar/sigma_bar) ./ ...
				factor.^((nu+1)/2);
		else
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
		end
		
	otherwise
		error('Invalid model : %s', model);
end


end

