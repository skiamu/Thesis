
function Statistics = PortfolioStatistics(Returns,freq)
%PortfolioStatistic given a multivariate returns times-eries computes
%several statistics
%
% INPUT:
%    Returns = Returns = matrix of asset class returns in the desired frequency
%    freq = desired return frequency, freq must be one of the following
%             string {'d','wk','m','q','s','y'} [string]
% OUTPUT:
%    Statistics.ExpReturnsAnn = expected returns annualized
%    Statistics.VolatilityAnn = volatility annualized
%    Statistics.Skew = Skewness
%    Statistics.Kurt = Kurtosis
%    Statistics.VaR = monthly value-at-risk (alpha = 5%)
%    Statistics.MaxDrawdown = maximum drawdown
%    Statistics.MeanDrawdown = mean drawdown
%    Statistics.Corr = Correlationmatrix

[~, n] = size(Returns);
%% 1) introductory plot
% return distribution histogram
nbins = 15;
figure
hist(Returns,nbins); % check return's distribution
if n > 1
	legend('Money','Bond','Equity')
	title('asset class returns histogram')
	xlabel('weekly returns')
% 	saveas(gcf,'/home/andrea/Thesis/Latex/final/Images/ReturnsHist.jpeg');
   print('/home/andrea/Thesis/Latex/final/Images/ReturnsHist', '-dpng', '-r900');
end
% returns plot
figure
for i = 1 : n
	subplot(n,1,i)
	plot(Returns(:,i))
	title(['returns time series',num2str(i)])
end

%% 2) asset class statistics
% 2.1) drawdowns
CumReturns = cumprod(1+Returns)-1; % cumulative return
AD = cummax(CumReturns) - CumReturns; % drawdown
MaxAD = max(AD); % maximum drawdown
MeanAD = mean(AD); % mean drawdown
figure
area(AD(:,3))
% hold on
% plot(Returns(:,3))

% 2.2) expected returns and volatility, annualized
switch freq
	case 'd'
		t = 252;
		t_VaR = 20;
	case 'wk'
		t = 52;
		t_VaR = 4;
	case 'm'
		t = 12;
		t_VaR = 1;
	case 'q'
		t = 4;
		t_VaR = 1 / 3;
	case 's'
		t = 2;
		t_VaR = 1 / 6;
	case 'y'
		t = 1;
		t_VaR = 1 / 12;
end
ExpReturnsAnn = mean(Returns) * t;
VolatilityAnn = std(Returns) * sqrt(t);
Median = median(Returns) * t;
Sharpe = (ExpReturnsAnn - ExpReturnsAnn(1)) ./ VolatilityAnn;
% 2.2) ex-post V@R, historical simulation
VaR = quantile(-Returns,0.95) * t_VaR; % 95 percent quantile loss distribution


% 2.4) Skweness and Kurtosis
Skew = skewness(Returns);
Kurt = kurtosis(Returns);

% 2.5) correlations
if n> 1
	Corr = corr(Returns);
end

%% 3) normality test
HZmvntest(Returns); % multivariate normality test

%% 4) print results
disp('%%%%%%%%%%%%%  Sample Statistics  %%%%%%%%%%%%%')
disp(['Means (ann) : ',num2str(ExpReturnsAnn)])
disp(['StDevs (ann) : ',num2str(VolatilityAnn)])
disp(['Median (ann) : ',num2str(Median)])
disp(['Skeweness : ',num2str(Skew)])
disp(['Kurtosis : ',num2str(Kurt)])
if n>1
	disp('Correlation matrix: ')
	Corr
end
disp(['Monthly VaR_0.95: ', num2str(VaR)])
disp(['Max Drawdown: ', num2str(MaxAD)])
disp(['Mean Drawdown: ', num2str(MeanAD)])
disp(['Sharpe Ratio: ', num2str(Sharpe)])

%% 5) return results
Statistics.ExpReturnsAnn = ExpReturnsAnn;
Statistics.VolatilityAnn = VolatilityAnn;
Statistics.Skew = Skew;
Statistics.Kurt = Kurt;
Statistics.VaR = VaR;
Statistics.MaxDrawdown = MaxAD;
Statistics.MeanDrawdown = MeanAD;
Statistics.Corr = Corr;

end % PortfolioStatistics




