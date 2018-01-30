
function Statistics = PortfolioStatistics(Returns,freq,policy,r,N)
%PortfolioStatistic computes several statistics on return data obtained by
%a given portfolio strategy 
%
% INPUT:
%    Returns = matrix of asset class returns in the desired frequency
%    freq = desired return frequency, freq must be one of the following
%             string {'d','wk','m','q','s','y'} [string]
%    policy = {'ODAA','CPPI','ConstantMix'}
%    r = annualized cash return (for the Sharpe index computation)
%    N = number of portfolio rebalancings
% OUTPUT:
%    Statistics.ExpReturnsAnn = expected returns annualized
%    Statistics.VolatilityAnn = volatility annualized
%    Statistics.Skew = Skewness
%    Statistics.Kurt = Kurtosis
%    Statistics.VaR = monthly value-at-risk (alpha = 5%)
%    Statistics.MaxDrawdown = maximum drawdown
%    Statistics.MeanDrawdown = mean drawdown
%    Statistics.Corr = Correlationmatrix


%% 1) asset class statistics
% 1.1) drawdowns
CumReturns = cumprod(1+Returns)-1; % cumulative return
AD = cummax(CumReturns) - CumReturns; % drawdown
MaxAD = max(AD); % maximum drawdown
MeanAD = mean(AD); % mean drawdown

% 1.2) expected returns and volatility, annualized
switch freq
	case 'd'
		t = N / 252;
		t_VaR = 20;
	case 'wk'
		t = N / 52;
		t_VaR = 4;
	case 'm'
		t = N / 12;
		t_VaR = 1;
	case 'q'
		t = N / 4;
		t_VaR = 1 / 3;
	case 's'
		t = N / 2;
		t_VaR = 1 / 6;
	case 'y'
		t = N / 1;
		t_VaR = 1 / 12;
end
InvestmentReturns = CumReturns(end,:); % the path are along the coluns
ExpReturnsAnn = (1 + mean(InvestmentReturns)).^(1/t) - 1;
VolatilityAnn = std(InvestmentReturns) * sqrt(1/t);
Median = (1 + median(InvestmentReturns)).^(1/t) - 1;
Sharpe = (ExpReturnsAnn - r) ./ VolatilityAnn;

% 1.3) ex-post V@R, historical simulation
VaR = (1 + quantile(-Returns,0.95)).^ t_VaR - 1; % 95 percent quantile loss distribution

% 1.4) Skweness and Kurtosis
Skew = skewness(Returns);
Kurt = kurtosis(Returns);

%% 2) plot results

ksdensity(InvestmentReturns);
hold on
% title([policy,' empirical investment return density'])

%% 5) return results
Statistics.ExpReturnsAnn = ExpReturnsAnn;
Statistics.VolatilityAnn = VolatilityAnn;
Statistics.Median = Median;
Statistics.Skew = Skew;
Statistics.Kurt = Kurt;
Statistics.VaR = VaR;
Statistics.MaxDrawdown = MaxAD;
Statistics.MeanDrawdown = MeanAD;
Statistics.Sharpe = Sharpe;
end % PortfolioStatistics




