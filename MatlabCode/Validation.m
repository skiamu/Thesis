function [p_star,Statistics] = Validation(w,X,U,Nsim,M,Nstep,freq,policy,r)
%Validation is a function for model validation. By giving in inputs asset
%allocation maps,target sets and simulated asset class returns this function
%builds up the paths and computes the probability of reaching the target set.
%   INPUT:
%      w = simulated asset class returns
%      X = discretized target sets for each time step (0,..,Nstep) [cell array]
%      U = asset allocation maps for each time step (0,..,Nstep-1) [cell array]
%      Nsim = number of MC simulation
%      M = asset allocation dimension
%      Nstep = number of time steps
%      freq =
%      policy = {'ODAA','CPPI','ConstantMix'}
%      r = annualized cash return (for the Sharpe index computation)
%   OUTPUT:
%      p_star = probability of reaching the target sets

%% Portfolio paths
xk = zeros([Nsim Nstep+1]);
xk(:,1) = X{1}; % initial porftolio value
u0 = U{1}'; % initial asset allocation
xk(:,2) = xk(:,1) .* (1 + w(:,:,1) * u0); % portfolio value at time k = 1
for k = 2 : Nstep
	Uinterp = zeros([Nsim M]); % maps matrix computed on portfolio realization
	idx1 = xk(:,k) > X{k}(end); % indexes ptf realizations greater than target set upper bound
	idx2 = xk(:,k) < X{k}(1); % indexes ptf realizations smaller than target set lower bound
	idx3 = ~(idx1 | idx2); % indexes ptf realization inside target sets
	if ~isempty(Uinterp(idx1,:))
		Uinterp(idx1,:) = repmat(U{k}(end,:), [sum(idx1) 1]);
	end
	if ~isempty(Uinterp(idx2,:))
		Uinterp(idx2,:) = repmat(U{k}(1,:), [sum(idx2) 1]);
	end
	[m,~]=size(U{k});
	if m ==1 % in the constant-mix case there's no need to interpolate
		xk(:,k+1) = xk(:,k) .* (1 + sum(U{k} .* w(:,:,k),2));
	else
		Uinterp(idx3,:) = interp1(X{k},U{k},xk(idx3,k)); % interpolate input maps on portfolio realizations
		xk(:,k+1) = xk(:,k) .* (1 + sum(Uinterp .* w(:,:,k),2)); % update portfolio value
	end
end
p_star = sum(xk(:,end) > X{end}(1)) / Nsim;

%% Results
% compure portfolio returns
Returns = xk(:,2:end) ./ xk(:,1:end-1) - 1;
Statistics = PortfolioStatistics(Returns',freq,policy,r,Nstep);
% print results
disp(['%%%%%%%%%%%%%  Sample Portfolio Return Statistics ',policy,' strategy  %%%%%%%%%%%%%'])
disp(['Means (ann) : ',num2str(mean(Statistics.ExpReturnsAnn))])
disp(['StDevs (ann) : ',num2str(mean(Statistics.VolatilityAnn))])
disp(['Median (ann) : ',num2str(mean(Statistics.Median))])
disp(['Skeweness : ',num2str(mean(Statistics.Skew))])
disp(['Kurtosis : ',num2str(mean(Statistics.Kurt))])
disp(['Monthly VaR_0.95: ', num2str(mean(Statistics.VaR))])
disp(['Max Drawdown: ', num2str(max(Statistics.MaxDrawdown))])
disp(['Mean Drawdown: ', num2str(mean(Statistics.MeanDrawdown))])
disp(['Sharpe Ratio: ', num2str(mean(Statistics.Sharpe))]) % da sistemare

end % validation

