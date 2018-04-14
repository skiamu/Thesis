function[U,Floor,Cushion] = CPPI(u0,X,r,m,N,param,model,VaR,alpha)
%CPPI is a function for implementing the CPPI strategy. It computed the
%allocation maps according to the CPPI policy for different realization of
%the portfolio value
%   INPUT:
%      u0 = initial portfolio allocation [column vector]
%      X = cell array of discretized target sets
%      r = expected cash return (in the right frequency)
%      m = CPPI multiplier 
%      N = number of time steps
%      param = cell array or struct of parameters
%      model =
%      VaR =
%      alpha =
%   OUTPUT:
%      U = asset allocation maps [cell array]
%      Floor = vector of portfolio floors, one at each time step
%      Cushion = cell array of portfolio cushions
%   REMARKS:

%% CPPI synthesis
M = length(u0); % asset class dimension
x0 = X{1}; % initial portfolio value
idxRiskyBasket = [2 3];
A = zeros([1 M]);
A(idxRiskyBasket) = 1;
beta = .9; % percentage guaranteed capital
Floor = zeros([N 1]); 
% Floor(1) = x0 * (1 - A * u0 / m); % floor that guarantees exposure E0
Floor(1) = beta * x0 / (1 + r)^N;
U = cell([N 1]); U{1} = u0';
Cushion = cell([N 1]); Cushion{1} = max(x0 - Floor(1),0);
Aeq = ones([1 M]); beq = 1;
ub = ones([M 1]); lb = zeros([M 1]);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
for k = 2 : N
	k
	Floor(k) = Floor(k-1) * (1 + r); % floor grows deterministically at cash rate
	dimXk = length(X{k});
	Uk = zeros([dimXk M]);
	for j = 1 : dimXk
		Cushion{k}(j) = max(X{k}(j) - Floor(k),0); % define portoflio cushion
		[Uk(j,:), ~] = fmincon(@(u)-objfun(u,A),u0,[],[],Aeq,beq,...
			lb,ub,@(u)confuneq(u,param,model,VaR,alpha,X{k}(j),Cushion{k}(j),m,A),options);
		u0 = Uk(j,:)';
	end
	U{k} = Uk;
end

end % CPPI

function f = objfun(u,A)
% the goal is to invest as much as possible in the risk asset (Equity asset
% class)
f = A * u;
end % objfun


function [c, ceq] = confuneq(u,param,model,VaR,alpha,x,Cushion,m,A)
%confuneq states the max problem's constraints
%   INPUT:
%      u = asset allocation vector
%      param = density parameters
%   OUTPUT:
%      c = inequality constraints
%      ceq = equality constraints
switch model
	case 'Gaussian'
		S = param.S; mu = param.mu;
		ceq = []; 
		mu_p = -u' * mu; sigma_p = sqrt(u' * S * u);
		c = [-VaR + (mu_p + norminv(1-alpha) * sigma_p);  % variance constraint
			A*u*x - m*Cushion];
	case 'Mixture'
		% 1) compute covariance matrix of asset return w(k+1)
		W = 0;
		for i = 1 : length(param)
			W = W + param(i).lambda * param(i).S;
			for j = 1 : i-1
				W = W + param(i).lambda * param(j).lambda * ...
					(param(i).mu - param(j).mu) * (param(i).mu - param(j).mu)';
			end
		end
		% 2) setup the constrints
		% By using this formula we suppose gaussianity, mu = 0 and
		% scaling rule, see how it can be improved
		sigma_max = VaR / norminv(1-alpha);
		ceq = []; 
		c = [u' * W * u - sigma_max^2; % variance constraint
			A*u*x - m*Cushion];  % exposure bounded by m*Cuschion
	case 'GH'
		% scrivere vincolo var in forma analitica
		lambda = param.lambda; Chi = param.Chi; Psi = param.Psi;
		Sigma = param.sigma; gamma = param.gamma;
		% 1) compute Cov[w(k+1)] = E[w(k+1)]*Sigma + Var[w(k+1)]*gamma*gamma'
		wMean = (Chi/Psi) * besselk(lambda+1,sqrt(Chi*Psi)) / besselk(lambda,sqrt(Chi*Psi));
		wMoment2 = (Chi/Psi)^2 * besselk(lambda+2,sqrt(Chi*Psi)) / besselk(lambda,sqrt(Chi*Psi));
		wVar = wMoment2 - wMean^2;
		wCov = wMean * Sigma + wVar * (gamma * gamma');
		W = wCov;
		% 2) setup the constrints
		% By using this formula we suppose gaussianity, mu = 0 and
		% scaling rule, see how it can be improved
		sigma_max = VaR / norminv(1-alpha);
		ceq = []; 
		c = [u' * W * u - sigma_max^2; % variance constraint
			A*u*x - m*Cushion];  % exposure bounded by m*Cuschion
	otherwise
		error('invalid model %s',model);
end


end







