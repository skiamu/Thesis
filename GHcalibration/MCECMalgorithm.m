function [theta,LogL,exitFlag,numIter] = MCECMalgorithm(toll,maxiter,X,GHmodel)
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
[mu,Sigma,gamma,lambda,alpha_bar] = InitialConditions(X,d);
X_bar = mean(X);
theta = cell([maxiter 1]);
theta{1} = {lambda,alpha_bar,mu,Sigma,gamma};
Q2 = zeros([maxiter 1]);
for k = 2 : maxiter
	k
	% 2) Compute delta, eta
	[delta,eta] = ComputeWeight(lambda,alpha_bar,d,X,mu,Sigma,gamma);
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
	[delta,eta,csi] = ComputeWeight(lambda,alpha_bar,d,X,mu,Sigma,gamma);
	
	% 6) maximize Q2(lambda,Chi,Psi)
	options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
	x0 = [lambda;alpha_bar];
	[A,b,Aeq,beq] = confuneq(GHmodel,d);
	[x_star,f_star] = fmincon(@(x)objfun(x,delta,eta,csi,GHmodel), ...
		x0,A,b,Aeq,beq,[],[],[],options);
	Q2(k) = -f_star;
	lambda = x_star(1); alpha_bar = x_star(2); 
	% go the to the first parametrization
	Psi = alpha_bar * besselk(lambda+1,alpha_bar) / besselk(lambda,alpha_bar);
   Chi = alpha_bar^2 / Psi;
	
	% 7) check conditions
	theta{k} = {lambda,alpha_bar,mu,Sigma,gamma};
	if checkCondition(toll,theta,k,Q2)
		exitFlag = 'condition';
		numIter = k;
		LogL = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma);
		theta{k}{end+1} = Chi; theta{k}{end+1} = Psi;
		break
	end
	numIter = k;
	LogL = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma);
	theta{k}{end+1} = Chi; theta{k}{end+1} = Psi;
end


end % MCECMalgorithm

function [delta,eta,csi] = ComputeWeight(lambda,alpha_bar,d,X,mu,Sigma,gamma)
%ComputeWeight is a function for computing the weigths delta, eta and csi.
%See formulas in the theory for analitic expression
%   INPUT:
%   OUTPUT:
%

% recover Chi and Psi from alpha_bar
Psi = alpha_bar * besselk(lambda+1,alpha_bar) / besselk(lambda,alpha_bar);
Chi = alpha_bar^2 / Psi;
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


function f = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,Gamma)
%LogLikelihood computed the loglikelihood function at maximum
%   INPUT:
%   OUTPUT:

[N, d] = size(X);
f = 0;
invSigma = inv(Sigma);

c = sqrt(Chi*Psi)^(-lambda) * Psi^lambda * (Psi+Gamma'*invSigma*Gamma)^(0.5*d - lambda) / ...
	((2*pi)^(d/2) * sqrt(det(Sigma)) * besselk(lambda,sqrt(Chi*Psi)));
for i = 1 : N
	x = X(i,:)';
	factor = sqrt((Chi + (x - mu)'*invSigma*(x - mu)) * (Psi + Gamma'*invSigma*Gamma));
	density = c * besselk(lambda-d/2,factor) .* exp((x-mu)'*invSigma*Gamma) ./ ...
		(factor).^(d/2 - lambda);
	f = f + log(density);
end

end % density


function [mu,Sigma,gamma,lambda,alpha_bar] = InitialConditions(X,d)
X_bar = mean(X); mu = X_bar';
Sigma = cov(X);
lambda = 1;
alpha_bar = 1;
gamma = zeros([d 1]);
end % initial conditions


