function [theta,LogL,exitFlag,numIter] = MCECMalgorithm2(toll,maxiter,X,GHmodel)
%MCECMalgorithm implements a modified version of the EM algorithm for
%fitting a Generalized Hyperbolic Distribution
%   INPUT:
%      toll = stopping tolerance
%      maxiter = maximum number of iterations
%      X = returns data
%      GHmodel = 't', 'VG', 'NIG'
%   OUTPUT:
%      theta = cell array of parameters
%      LogL = log-likelihood at optimum
%      exitFlag = 'maxiter' or 'condition'
%      numIter = number of algorithm iterations
% REMARKS: scrivere la funzione obiettivo Q2 anke per gli altri modelli

exitFlag = 'maxiter';
% 1) set initial conditions
[N, d] = size(X);
[mu,Sigma,gamma,lambda,Chi,Psi] = InitialConditions(X,GHmodel,d);
X_bar = mean(X);
theta = cell([maxiter 1]);
theta{1} = {lambda,Chi,Psi,mu,Sigma,gamma};
Q2 = zeros([maxiter 1]);
for k = 2 : maxiter
	k
	% 2) Compute delta, eta
	[delta,eta] = ComputeWeight(lambda,Chi,Psi,d,X,mu,Sigma,gamma,GHmodel);
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
	Sigma = det(Sigma)^(1/d) * D / det(D)^(1/d); % take the root with lowest phase, complex
% 	Sigma = nthroot(det(Sigma),d) * D / nthroot(det(D),d);
	% 5) update the weights
	[delta,eta,csi] = ComputeWeight(lambda,Chi,Psi,d,X,mu,Sigma,gamma,GHmodel);
	
	% 6) maximize Q2(lambda,Chi,Psi)
	options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
	x0 = [lambda;Chi;Psi];
	[A,b,Aeq,beq] = confuneq(GHmodel);
	[x_star,f_star] = fmincon(@(x)objfun(x,delta,eta,csi,GHmodel), ...
		x0,A,b,Aeq,beq,[],[],[],options);
	Q2(k) = -f_star;
	lambda = x_star(1); Chi = x_star(2); Psi = x_star(3);
	
	% 7) check conditions
	theta{k} = {lambda,Chi,Psi,mu,Sigma,gamma};
	if checkCondition(toll,theta,k,Q2)
		exitFlag = 'condition';
		numIter = k;
		LogL = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma,GHmodel);
		break
	end
	numIter = k;
	LogL = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma,GHmodel);
end


end % MCECMalgorithm

function [delta,eta,csi] = ComputeWeight(lambda,Chi,Psi,d,X,mu,Sigma,gamma,GHmodel)
%ComputeWeight is a function for computing the weigths delta, eta and csi.
%See formulas in the theory for analitic expression
%   INPUT:
%   OUTPUT:
%
[N,~] = size(X);
% 1) compute lambda_bar, Psi_bar and Chi_bar according to the model
lambda_bar = lambda - 0.5 * d;
invSigma = inv(Sigma);
Psi_bar =  gamma' * invSigma * gamma;
Chi_bar = zeros([N 1]);
for i = 1 : N % check if Sigma is the right matrix
	Chi_bar(i) = (X(i,:)' - mu)' * invSigma * (X(i,:)' - mu);
end
switch GHmodel
	case 't'
		Chi_bar = Chi_bar + Chi;
	case 'VG'
		Psi_bar = Psi + Psi_bar;
	otherwise % NIG, general case
		Psi_bar = Psi_bar + Psi;
		Chi_bar = Chi_bar + Chi;
end
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


function f = objfun(x,delta,eta,csi,GHmodel)
%objfun is the objective function Q2 for the maximization problem
%   INPUT:
%      x = vector of parameters, [lambda,Chi,Psi]
%      delta =
%      eta =
%      csi =
%   OUTPUT:
%      f =
% REMARKS: different objective function for 't' and 'VG'
[N, ~] = size(delta);
switch GHmodel
	case 't'
		f = -x(1) * N * log(x(2) / 2) - N * log(gamma(-x(1))) + (x(1)-1) * sum(csi) - ...
			0.5 * x(2) * sum(delta);
	case 'VG'
		f = x(1) * N * log(x(3) / 2) - N * log(gamma(x(1))) + (x(1) - 1) * sum(csi) - ...
			0.5 * x(3) * sum(eta);
	otherwise
		f = (x(1)-1) * sum(csi) - 0.5 * x(2) * sum(delta) - 0.5 * x(3) * sum(eta) -...
			0.5 * N * x(1) * log(x(2)) + 0.5 * N * x(1) * log(x(3)) - ...
			N * log(2 * besselk(x(1),sqrt(x(2) * x(3))));
end

f = - f; % maximization problem
end % objfun


function [A,b,Aeq,beq] = confuneq(GHmodel)
%confuneq is a function for setting the constraint for the maximization
%problem
%   INPUT:
%   OUTPUT:

switch GHmodel
	case 'NIG' % lambda = 0.5
		A = [0 -1 0; 0 0 -1]; b = [eps; 0];
		Aeq = [1 0 0]; beq = -0.5;
	case 'VG'
		% Chi = 0, Psi > 0 lambda > 0
		A = [-1 0 0; 0 0 -1]; b = [eps; eps];
		Aeq = [0 1 0]; beq = 0;
	case 't'
		% lambda < 0, Psi = 0, Chi > 0
		Aeq = [0 0 1]; beq = 0;
		A = [0 -1 0; 1 0 0]; b = [eps; eps];
	otherwise
		Aeq = []; beq = [];
		A = [0 -1 0; 0 0 -1]; b = [eps; eps];
end

end % objfun

function condition = checkCondition(toll,theta,k,Q2)
% REMARK: set different tollerance
[N,~] = size(theta{k});
diff = zeros([N 1]);
for i = 1 : N
	diff(i) = norm(theta{k}{i}-theta{k-1}{i});
end
if all(diff) < toll && norm(Q2(k) - Q2(k-1)) < toll
	condition = true;
else
	condition = false;
end
end % checkCondition

function f = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,Gamma,GHmodel)
%LogLikelihood computed the loglikelihood function at maximum
%   INPUT:
%   OUTPUT:

[N, d] = size(X);
f = 0;
invSigma = inv(Sigma);
switch GHmodel
	case 't'
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
	case 'VG'
		c = 2*lambda^(lambda)*(2*lambda+Gamma'*invSigma*Gamma)^(d/2-lambda) / ...
			((2*pi)^(d/2)*sqrt(det(Sigma))*gamma(lambda));
		for i = 1 : N
			x = X(i,:)';
			factor = sqrt((x-mu)'*invSigma*(x-mu)*(2*lambda + Gamma'*invSigma*Gamma));
			density = c * besselk(lambda-d/2,factor) .* exp((x-mu)'*invSigma*Gamma) ./ ...
				factor.^(d/2 - lambda);
			f = f + log(density);
		end
	otherwise % NIG, general case
		c = sqrt(Chi*Psi)^(-lambda) * Psi^lambda * (Psi+Gamma'*invSigma*Gamma)^(0.5*d - lambda) / ...
			((2*pi)^(d/2) * sqrt(det(Sigma)) * besselk(lambda,sqrt(Chi*Psi)));
		for i = 1 : N
			x = X(i,:)';
			factor = sqrt((Chi + (x - mu)'*invSigma*(x - mu)) * (Psi + Gamma'*invSigma*Gamma));
			density = c * besselk(lambda-d/2,factor) .* exp((x-mu)'*invSigma*Gamma) ./ ...
				(factor).^(d/2 - lambda);
			f = f + log(density);
		end
end

end % density

function [mu,Sigma,gamma,lambda,Chi,Psi] = InitialConditions(X,GHmodel,d)
X_bar = mean(X); mu = X_bar';
Sigma = cov(X);
switch GHmodel
	case 't'
		gamma = ones([d 1])/1e5;
		lambda = -0.5; Chi = 0.7; Psi = 0.7;
	case 'VG'
		gamma = zeros([d 1]);
		lambda = 1; Chi = 0; Psi = 0.7;
	otherwise
		gamma = zeros([d 1]);
		lambda = 0.5; Chi = 0.7; Psi = 0.7;
end
end % initial conditions


