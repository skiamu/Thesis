function [ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha)
%DPalgorithm implements a Dynamic Programming algorithm to solve a
%stochastic reachability problem
%   INPUT:
%      N = number of time steps
%      M = dimension asset allocation (e.g. 3)
%      X = cell array of discretized target sets, to access the i-th
%          element use X{i}
%      param = density parameters (struct or cell array)
%      VaR = yearly value at risk
%      alpha = 
%   OUTPUT:
%      U = cell array of asset allocations
%   Remarks:
%      1) when the mixture is used, param is a cell array of the following
%         form {{mu1, Sigma1, lambda1}, ..., {mu_n, Sigma_n, lambda_n}}

if nargin < 8
	freq = 1; % weekly frequency
end
U = cell([N 1]); % cell array, in each position there'll be a matrix
J = cell([N+1 1]); % optimal value function
J{end} = ones([length(X{end}) 1]);  % indicator function of the last target set
% options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
for k = N : -1 : 1
	k
	u0 = ones([M 1])/M; % equally-weighted portfolio, initial guess, (maybe better past solution)
	dimXk = length(X{k});
	Uk = zeros([dimXk M]);
	Jk = zeros([dimXk 1]);
	for j = 1 : dimXk
		[Uk(j,:), Jk(j)] = fmincon(@(u)objfun(u,X{k}(j),X{k+1},J{k+1},param,model), ...
			u0,[],[],[],[],[],[],@(u)confuneq(u,param,model,VaR,alpha),options);
	end
	U{k} = Uk; J{k} = -Jk;
end
end %  DPalgorithm

function f = objfun(u,x,int_domain,J,param,model)
%objfun defines the problem's objective function that will be maximized
%   INPUT:
%      u = asset allocation vector
%      x = realization portfolio value [scalar]
%      int_domain = integration domain, k+1 target set discretization
%      J = k+1 optimal value function
%      param = density parameters
%   OUTPUT:
%      f = objective function

f = trapz(int_domain, J .* pf(int_domain,x,u,param,model));
f = -f; % solve the assocoated max problem

end % objfun

function [c, ceq] = confuneq(u,param,model,VaR,alpha)
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
		ceq = u' * ones(size(u)) - 1; % budget constraint
		mu_p = -u' * mu; sigma_p = u' * S * u;
		c = [-u;                                        % long-only constraint
			VaR - (mu_p + norminv(1-alpha) * sigma_p)];  % variance constraint
	case 'Mixture'
		n = length(param); 
		% 1) compute the cdf
		Phi = 0;
		for i = 1 : n
			mu = u' * param{i}{1};
			sigma = sqrt(u' * param{i}{2} * u);
			lambda = param{i}{3};
			Phi = Phi + lambda * normcdf((-VaR - mu) / sigma);
		end
		% 2) setup the constraints
		ceq = u' * ones(size(u)) - 1; % budget constraint
		c = [-u;            % long-only constraint
			Phi - alpha];    % variance constraint
	case 'Levy'
		
	otherwise
		error('invalid model %s',model);
end

end % confuneq
















