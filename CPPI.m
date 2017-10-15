function[U,Floor,Cushion] = CPPI(u0,X,r,m,N,param,model,VaR,alpha)
%CPPI is a function for implementing the CPPI strategy. It computed the
%allocation maps according to the CPPI policy for different realization of
%the portfolio value
%   INPUT:
%      u0 = initial portfolio allocation [column vector]
%      X = cell array of discretized target sets
%      r = expected cash return
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
idxRiskyBasket = 3;
Floor = zeros([N 1]);
Floor(1) = x0 * (1 - u0(idxRiskyBasket) / m); % floor that guarantees exposure E0
U = cell([N 1]); U{1} = u0';
Cushion = cell([N 1]); Cushion{1} = subplus(x0 - Floor(1));
options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
for k = 2 : N
	k
	Floor(k) = Floor(k-1) * (1 + r); % floor grows deterministically at cash rate
	dimXk = length(X{k});
	Uk = zeros([dimXk M]);
	for j = 1 : dimXk
		Cushion{k}(j) = subplus(X{k}(j) - Floor(k));
		[Uk(j,:), ~] = fmincon(@(u)objfun(u,idxRiskyBasket,M), ...
			u0,[],[],[],[],[],[],@(u)confuneq(u,param,model,VaR,alpha,X{k}(j),...
			Cushion{k}(j),m,idxRiskyBasket,M),options);
		u0 = Uk(j,:)';
	end
	U{k} = Uk;
end

end % CPPI

function f = objfun(u,indexRiskyBasket,M)
% the goal is to invest as much as possible in the risk asset (Equity asset
% class)
A = zeros([1 M]);
A(indexRiskyBasket) = 1;
f = A * u;
f = -f;
end % objfun


function [c, ceq] = confuneq(u,param,model,VaR,alpha,x,Cushion,m,indexRiskyBasket,M)
%confuneq states the max problem's constraints
%   INPUT:
%      u = asset allocation vector
%      param = density parameters
%   OUTPUT:
%      c = inequality constraints
%      ceq = equality constraints
VaRConstraintType = '2';
A = zeros([1 M]); A(indexRiskyBasket) = 1;
switch model
	case 'Gaussian'
		S = param.S; mu = param.mu;
		ceq = u' * ones(size(u)) - 1; % budget constraint
		mu_p = -u' * mu; sigma_p = sqrt(u' * S * u);
		c = [-u;                                        % long-only constraint
			-VaR + (mu_p + norminv(1-alpha) * sigma_p);  % variance constraint
			A*u*x - m*Cushion];  
	case 'Mixture'
		n = length(param);
		switch VaRConstraintType
			case '1' % analytic method
				
			case '2' % quadratic method
				% 1) compute covariance matrix of asset return w(k+1)
				W = zeros(length(u));
				for i = 1 : n
					W = W + param{i}{3} * param{i}{2};
					for j = 1 : i-1
						W = W + param{i}{3}*param{j}{3}*(param{i}{1}-param{j}{1})...
							*(param{i}{1}-param{j}{1})';
					end
				end
				% 2) setup the constrints
				% By using this formula we suppose gaussianity, mu = 0 and
				% scaling rule, see how it can be improved
				sigma_max = VaR / norminv(1-alpha);
				ceq = u' * ones(size(u)) - 1; % budget constraint
				c = [-u;            % long-only constraint
					u' * W * u - sigma_max^2; % variance constraint
					A*u*x - m*Cushion];  % exposure bounded by m*Cuschion
		end
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
		ceq = u' * ones(size(u)) - 1; % budget constraint
		c = [-u;            % long-only constraint
			u' * W * u - sigma_max^2; % variance constraint
			A*u*x - m*Cushion];  % exposure bounded by m*Cuschion
	otherwise
		error('invalid model %s',model);
end


end







