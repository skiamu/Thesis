function [p_star] = Validation(w,X,U,Nsim,M,Nstep)
%Validation is a function for model validation. By giving in inputs maps,
%target sets and simulated asset class returns this function builds up the
%paths and computes the probability of reaching the target set.
%   INPUT:
%      w = simulated asset class returns
%      X = discretized target sets for each time step (0,..,Nstep) [cell array]
%      U = asset allocation maps for each time step (0,..,Nstep-1) [cell array]
%      Nsim = number of MC simulation
%      M = asset allocation dimension
%      Nstep = number of time steps
%   OUTPUT:
%      p_star = probability of reaching the target sets

%% Portfolio paths
xk = X{1}; % initial porftolio value
u0 = U{1}'; % initial asset allocation
xk = xk * (1 + w(:,:,1) * u0); % portfolio value at time k = 1
for k = 2 : Nstep
	Uinterp = zeros([Nsim M]); % maps matrix computed on portfolio realization
	idx1 = xk > X{k}(end); % indexes ptf realizations greater than target set upper bound
	idx2 = xk < X{k}(1); % indexes ptf realizations smaller than target set lower bound
	idx3 = ~(idx1 | idx2); % indexes ptf realization inside target sets
	if ~isempty(Uinterp(idx1,:))
		Uinterp(idx1,:) = U{k}(end,:);
	end
	if ~isempty(Uinterp(idx2,:))
		Uinterp(idx2,:) = U{k}(1,:);
	end
	Uinterp(idx3,:) = interp1(X{k},U{k},xk(idx3)); % interpolate input maps on portfolio realizations
	xk = xk .* (1 + sum(Uinterp .* w(:,:,k),2)); % update portfolio value
end
p_star = length(xk(xk > X{end}(1))) / Nsim;

%% Results



end % validation

