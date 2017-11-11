function [param, CalibrationData] = MCECMalgorithm_t(toll,maxiter,X,GHmodel)
%MCECMalgorithm implements a modified version of the EM algorithm for
%fitting a Generalized Hyperbolic Distribution
%   INPUT:
%      toll = stopping tolerance
%      maxiter = maximum number of iterations
%      X = returns data
%      GHmodel = 't', 'VG', 'NIG'
%   OUTPUT:
%      param = 
%      CalibrationData = 
% REMARKS: scrivere la funzione obiettivo Q2 anke per gli altri modelli

exitFlag = 'maxiter';
% 1) set initial conditions
[N, d] = size(X);
[mu,Sigma,gamma,lambda] = InitialConditions(X,d);
X_bar = mean(X);
theta = cell([maxiter 1]);
theta{1} = {lambda,mu,Sigma,gamma};
Q2 = zeros([maxiter 1]);
for k = 2 : maxiter
	k
	% 2) Compute delta, eta
	[delta,eta] = ComputeWeight(lambda,d,X,mu,Sigma,gamma);
	deltaMean = sum(delta) / N;
	etaMean = sum(eta) / N;
	
	% 3) update gamma
	gamma = (sum(delta * ones([1 d]) .* (repmat(X_bar,[N 1]) - X))' / N) / ...
		(deltaMean * etaMean - 1);
	
	% 4) update mu and Sigma
	mu = (sum(delta * ones([1 d]) .* X)' / N - gamma) / deltaMean;
	D = 0;
	for i = 1 : N
		D = D + delta(i) * (X(i,:)' - mu) * (X(i,:)' - mu)';
	end
	D = D/N - etaMean * (gamma * gamma');
	Sigma = D;
	% 5) update the weights with new mu,sigma and gamma
	[delta,eta,csi] = ComputeWeight(lambda,d,X,mu,Sigma,gamma);
	
	% 6) maximize Q2(lambda,Chi,Psi)
	options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
	x0 = lambda;
	[A,b,Aeq,beq] = confuneq(GHmodel,d);
	[x_star,f_star] = fmincon(@(x)objfun(x,delta,eta,csi,GHmodel), ...
		x0,A,b,Aeq,beq,[],[],[],options);
	Q2(k) = -f_star;
	lambda = x_star;
	% go the to the first parametrization
	nu = -2 * lambda;
	
	% 7) check conditions
	theta{k} = {lambda,mu,Sigma,gamma};
	if checkCondition(toll,theta,k,Q2)
		exitFlag = 'condition';
		numIter = k;
		LogL = LogLikelihood(X,lambda,mu,Sigma,gamma);
		theta{k}{end+1} = nu;
		break
	end
	numIter = k;
	LogL = LogLikelihood(X,lambda,mu,Sigma,gamma);
	theta{k}{end+1} = nu;
end
% 8) build the param struct
param.lambda = theta{numIter}{1};
param.mu = theta{numIter}{2};
param.sigma = theta{numIter}{3};
param.gamma = theta{numIter}{4};
param.nu = theta{numIter}{5};
nParam =  1 + d + d*(d+1) / 2 + d; % t
AIC = -2 * LogL + 2 * nParam;

% 9) return calibration information
CalibrationData.LogL = LogL;
CalibrationData.ExitFlag = exitFlag;
CalibrationData.numIter = numIter;
CalibrationData.AIC = AIC;

end % MCECMalgorithm_t

function [delta,eta,csi] = ComputeWeight(lambda,d,X,mu,Sigma,gamma)
%ComputeWeight is a function for computing the weigths delta, eta and csi.
%See formulas in the theory for analitic expression
%   INPUT:
%   OUTPUT:
%

% recover Chi and Psi from alpha_bar
nu = -2 * lambda;
Chi = nu - 2;
Psi = 0;
[N,~] = size(X);
% 1) compute lambda_bar, Psi_bar and Chi_bar according to the model NIG,
%    generelized hyperbolic and hyperbolic
lambda_bar = lambda - 0.5 * d;
invSigma = inv(Sigma);
Psi_bar =  gamma' * invSigma * gamma + Psi;
Chi_bar = zeros([N 1]);
for i = 1 : N % check if Sigma is the right matrix
	Chi_bar(i) = (X(i,:)' - mu)' * invSigma * (X(i,:)' - mu);
end
Chi_bar = Chi_bar + Chi;

% 2) compute weights
delta = (Chi_bar / Psi_bar).^(-0.5) .* besselk(lambda_bar-1,sqrt(Chi_bar * Psi_bar)) ./ ...
	besselk(lambda_bar,sqrt(Chi_bar * Psi_bar));
eta = (Chi_bar / Psi_bar).^(0.5) .* besselk(lambda_bar+1,sqrt(Chi_bar * Psi_bar)) ./ ...
	besselk(lambda_bar,sqrt(Chi_bar * Psi_bar));
if nargout > 2
	h = 0.001;
	csi = ((Chi_bar / Psi_bar).^(h/2) .* besselk(lambda_bar+h,sqrt(Chi_bar * Psi_bar)) ./ ...
		besselk(lambda_bar,sqrt(Chi_bar * Psi_bar)) - 1) / h;
end

end %ComputeWeight

function f = LogLikelihood(X,lambda,mu,Sigma,Gamma)
%LogLikelihood computed the loglikelihood function at maximum
%   INPUT:
%   OUTPUT:

[N, d] = size(X);
f = 0;
invSigma = inv(Sigma);

nu = -2*lambda; % degrees of freedom
c = (nu-2)^(nu/2) * (Gamma'*invSigma*Gamma)^((nu+d)/2) / ...
	((2*pi)^(d/2)*sqrt(det(Sigma))*gamma(nu/2)*2^(nu/2-1));
for i = 1 : N
	x = X(i,:)';
	factor = sqrt((nu - 2 + (x-mu)'*invSigma*(x-mu)) * Gamma'*invSigma*Gamma);
	density = c * besselk((nu+d)/2,factor) .* exp((x-mu)'*invSigma*Gamma) ./ ...
		factor.^((nu + d)/2);
	f = f + log(density);
end
end % density


function condition = checkCondition(toll,theta,k,Q2)
% REMARK: set different tollerance
[N,~] = size(theta{k});
diff = zeros([N 1]);
for i = 1 : N
	diff(i) = norm(theta{k}{i}-theta{k-1}{i});
end
if (all(diff) < toll && norm(Q2(k) - Q2(k-1)) < toll) 
	condition = true;
else
	condition = false;
end
end % checkCondition


function [mu,Sigma,gamma,lambda] = InitialConditions(X,d)
X_bar = mean(X); mu = X_bar';
Sigma = cov(X);
lambda = -2;
gamma = ones([d 1])/1e5;
end % initial conditions


