function [theta,LogL,exitFlag,numIter] = MCECMalgorithm(toll,maxiter,X,GHmodel)
%MCECMalgorithm implements a modovied version of the EM algorithm for
%fitting a Generalized Hyperbolic Distribution
%   INPUT:
%   OUTPUT:
% REMARKS: scrivere la funzione obiettivo Q2 anke per gli altri modelli
exitFlag = 'maxiter';
% 1) set initial conditions
[N, d] = size(X);
X_bar = mean(X);
mu = X_bar'; Sigma = cov(X); gamma = zeros([d 1]);
lambda = 0.5; Chi = 7; Psi = 7;
theta = cell([maxiter 1]);
theta{1} = {lambda,Chi,Psi,mu,Sigma,gamma};
Q2 = zeros([maxiter 1]);
for k = 2 : maxiter
	k
	% 2) Compute delta, eta
	[delta,eta] = ComputeWeight(lambda,Chi,Psi,d,X,mu,Sigma,gamma);
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
	% 	Sigma = det(Sigma)^(1/d) * D / det(D)^(1/d); % take the root with lowest
	% 	phase, complex
	Sigma = nthroot(det(Sigma),d) * D / nthroot(det(D),d);
	% 5) update the weights
	[delta,eta,csi] = ComputeWeight(lambda,Chi,Psi,d,X,mu,Sigma,gamma);
	
	% 6) maximize Q2(lambda,Chi,Psi)
	options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
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
		LogL = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma);
		break
	end
	numIter = k;
	LogL = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma);
end


end % MCECMalgorithm

function [delta,eta,csi] = ComputeWeight(lambda,Chi,Psi,d,X,mu,Sigma,gamma)
%ComputeWeight is a function for computing the weigths delta, eta and csi.
%See formulas in the theory for analitic expression
%   INPUT:
%   OUTPUT:
%
[N,~] = size(X);
% 1) compute lambda_bar, Psi_bar and Chi_bar
lambda_bar = lambda - 0.5 * d;
Psi_bar = Psi + gamma' * inv(Sigma) * gamma;
Chi_bar = zeros([N 1]);
for i = 1 : N % check if Sigma is the right matrix
	Chi_bar(i) = (X(i,:)' - mu)' * inv(Sigma) * (X(i,:)' - mu);
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
	case 'VG'
		f = N * x(1) * log(0.5 * x(3)) - N * log(gamma(x(1))) + (x(1) - 1) * sum(csi) - ...
			0.5 * sum(eta);
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
	case 'NIG'
		A = [0 -1 0; 0 0 -1]; b = [eps; 0];
		Aeq = [1 0 0]; beq = -0.5;
	case 'VG'
		% Chi = 0, Psi > 0 lambda > 0
		A = [0 1 0]; b = 0;
		Aeq = [0 0 -1;-1 0 0]; beq = [eps; eps];
	case 't'
		A = [0 0 1]; b = 0;
		Aeq = [0 -1 0]; beq = 0;
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

function f = LogLikelihood(X,lambda,Chi,Psi,mu,Sigma,gamma)
%LogLikelihood computed the loglikelihood function at maximum
%   INPUT:
%   OUTPUT:

[N, d] = size(X);
f = 0;
c = sqrt(Chi*Psi)^(-lambda) * Psi^lambda * (Psi+gamma'*inv(Sigma)*gamma)^(0.5*d - lambda) / ...
	((2*pi)^(d/2) * sqrt(det(Sigma)) * besselk(lambda,sqrt(Chi*Psi)));
for i = 1 : N
	x = X(i,:)';
	factor = sqrt((Chi + (x - mu)'*inv(Sigma)*(x - mu)) * (Psi + gamma'*inv(Sigma)*gamma));
	density = c * besselk(lambda-d/2,factor) .* exp((x-mu)'*inv(Sigma)*gamma) ./ ...
		(factor)^(d/2 - lambda);
	f = f + log(density);
end
end % density




