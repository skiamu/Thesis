function [ U, J] = DPalgorithm(N,M,X,param)
%DPalgorithm implements a Dynamic Programming algorithm to solve a
%stochastic reachability problem
%   INPUT:
%      N = number of time steps
%      M = dimension asset allocation (e.g. 3)
%      X = cell array of discretized target sets, to access the i-th
%          element use X{i}
%      param = vector of density parameters
%   OUTPUT:
%      U = cell array of asset allocations

U = cell([N 1]); % cell array, in each position there'll be a matrix
J = cell([N+1 1]); 
J{end} = ones([length(X{end}) 1]);  % indicator function of the last target set
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
for k = N : -1 : 1
	u0 = ones([M 1])/M; % equally-weighted portfolio, initial guess, (maybe better past solution)
	dimXk = length(X{k});
	Uk = zeros([dimXk M]);
	Jk = zeros([dimXk 1]);
	for j = 1 : dimXk 
		[Uk(j,:), Jk(j)] = fmincon(@(u)objfun(u,X{k}(j),X{k+1},J{k+1},param), ...
			u0,[],[],[],[],[],[],@(u)confuneq(u,param),options);
	end
	U{k} = Uk; J{k} = -Jk;
end
end %  DPalgorithm

function f = objfun(u,x,int_domain,J,param)
%objfun defines the problem's objective function that will be maximized
%   INPUT:
%      u = asset allocation vector
%      x = realization portfolio value [scalar]
%      int_domain = integration domain, k+1 target set discretization
%      J = k+1 optimal value function
%      param = density parameters
%   OUTPUT:
%      f = 
f = trapz(int_domain, J .* pf(int_domain,x,u,param));
f = -f; % solve the assocoated max problem
end % objfun

function [c, ceq] = confuneq(u,param)
%confuneq states the max problem's constraints
%   INPUT:
%      u = 
%      para = 
%   OUTPUT:
%      c = 
%      ceq = 
% S = param.S;
% rb = 0.0017;
% ceq = [u' * ones(size(u)) - 1; % budget constraint
% 	u' * S * u - rb]; % variance constraint
ceq = u' * ones(size(u)) - 1;
c = -u; % long-only constraint
end
















