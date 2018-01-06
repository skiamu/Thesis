function [p_star,Times] = ValidationDES(X,U,Nsim,Nstep,param,J_jump)
%ValidationDES is a function for the validation of the model in a
%event-driven approach
%   INPUT:
%      X = discretized target sets for each time step (0,..,Nstep) [cell array]
%      U = asset allocation maps for each time step (0,..,Nstep-1) [cell array]
%      Nsim = number of MC simulation
%      Nstep = number of time steps
%      lambda = speed parameter
%      p = trend parameter
%      r = cash yearly interest rate (cont compounded)
%      J_jump = jump size
%   OUTPUT:
%      p_star = probability of reaching the target sets

% parameters extraction
p = param.p; lambda = param.lambda;r = param.r;
%% 1) simulation random variables
% 1.1) simulation exponential r.v.
tau = exprnd(1 / lambda, [Nsim Nstep]);
% 1.2) simulation bernoulli r.v.
deltaN = -ones([Nsim Nstep]);
deltaN(rand([Nsim Nstep]) < p) = 1;

%% 2) validation
xk = X{1}; % initial wealth
u0 = U{1}; % initial asset allocation
xk = xk * (exp(r * tau(:,1)) + u0 * J_jump * deltaN(:,1));
for k = 2 : Nstep
	Uinterp = zeros([Nsim 1]);
	idx1 = xk > X{k}(end); % indexes ptf realizations greater than target set upper bound
	idx2 = xk < X{k}(1); % indexes ptf realizations smaller than target set lower bound
	idx3 = ~(idx1 | idx2); % indexes ptf realization inside target sets
	if ~isempty(Uinterp(idx1))
		Uinterp(idx1) = U{k}(end);
	end
	if ~isempty(Uinterp(idx2))
		Uinterp(idx2) = U{k}(1);
	end
	Uinterp(idx3) = interp1(X{k},U{k},xk(idx3));
	xk = xk .* (exp(r * tau(:,k)) + Uinterp .* J_jump .* deltaN(:,k));
end
p_star = length(xk(xk > X{end}(1))) / Nsim;
Times = sum(tau,2); % time to reach the target in years

end % ValidationDES

