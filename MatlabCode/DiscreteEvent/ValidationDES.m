function [p_star,Times] = ValidationDES(X,U,Nsim,Nstep,param,J_jump,model,theta,n)
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
%      model = {basic,ext1,ext2}
%   OUTPUT:
%      p_star = probability of reaching the target sets
rng default
r = param.r;
%% 1) simulation random variables
[Binomial, tau] = SimulationED(param,Nsim, Nstep,J_jump,model);

%% 2) validation
xk = zeros([Nsim Nstep+1]);
xk(:,1) = X{1}; % initial porftolio value
u0 = U{1}; % initial asset allocation
xk(:,2) = xk(:,1) .* (exp(r * tau(:,1)) + u0 * J_jump * Binomial(:,1));
for k = 2 : Nstep
	Uinterp = zeros([Nsim 1]);
	idx1 = xk(:,k) > X{k}(end); % indexes ptf realizations greater than target set upper bound
	idx2 = xk(:,k) < X{k}(1); % indexes ptf realizations smaller than target set lower bound
	idx3 = ~(idx1 | idx2); % indexes ptf realization inside target sets
	if ~isempty(Uinterp(idx1))
		Uinterp(idx1) = U{k}(end);
	end
	if ~isempty(Uinterp(idx2))
		Uinterp(idx2) = U{k}(1);
	end
	Uinterp(idx3) = interp1(X{k},U{k},xk(idx3,k));
	xk(:,k+1) = xk(:,k) .* (exp(r * tau(:,k)) + Uinterp .* J_jump .* Binomial(:,k));
end
p_star = sum(xk(:,end) > X{end}(1)) / Nsim;
Times = sum(tau,2); % time to reach the target in years

%% 3) Compute statistics
InvestmentReturns = xk(:,end) ./ xk(:,1) - 1;
InvestmentReturnsAnn = (1+InvestmentReturns).^(1./Times)-1;
ExpReturnsAnn = mean(InvestmentReturnsAnn);
VolatilityAnn = std(InvestmentReturnsAnn);
Median = median(InvestmentReturnsAnn);
Sharpe = (ExpReturnsAnn - (exp(r)-1)) ./ VolatilityAnn;
Skew = skewness(InvestmentReturnsAnn);
Kurt = kurtosis(InvestmentReturnsAnn);
% VaR_monthly = quantile(-InvestmentReturnsAnn,0.99);
[VaR,ES] = hHistoricalVaRES(InvestmentReturnsAnn,.99);
% R = xk(:,2:end)./xk(:,1:end-1)-1;
% CumR = cumprod(1+R,2)-1;
% AD = cummax(CumR,2)-CumR;
% MaxAD = max(AD,[],2);
%% 4) print results

% print results
disp('%%%%%%%%%%%%%  Sample Portfolio Return Statistics  %%%%%%%%%%%%%')
disp(['Mean Return (ann) : ',num2str(ExpReturnsAnn)])
disp(['Volatility (ann) : ',num2str(VolatilityAnn)])
disp(['Median Return (ann) : ',num2str(Median)])
disp(['Skeweness : ',num2str(Skew)])
disp(['Kurtosis : ',num2str(Kurt)])
disp(['Sharpe Ratio: ', num2str(Sharpe)]) % da sistemare
disp(['Average investment horizon: ',num2str(mean(Times))]); 
disp(['yearly VaR: ',num2str(VaR)]);
disp(['yearly ES: ',num2str(ES)]);
% disp(['Max DrawDown: ',num2str(mean(MaxAD))]);

figure
ksdensity(InvestmentReturnsAnn);
% yL = get(gca,'YLim');
% line([theta theta],yL,'Color','r');
title('investment return empirical density');
xlim([-0.2 0.2])
xlabel('Investent return');
print(strcat('/home/andrea/Thesis/Latex/final/Images/Density',model),'-dpng', '-r900');
figure
hist(InvestmentReturnsAnn,200)
title('Histogram investment return');
xlabel('Investent returns');
xlim([-0.2 0.2])
% yL = get(gca,'YLim');
% line([theta theta],yL,'Color','r');
print(strcat('/home/andrea/Thesis/Latex/final/Images/DensityHist',model),'-dpng', '-r900');

end % ValidationDES

