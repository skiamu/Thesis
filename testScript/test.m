clc;close all; clear variables;
freq = 'wk';
tickerMM = 'SHV';
tickerBM = 'BGBRX';
tickerEM = '^GSPC';
InitialDate = '23012000';
EndDate = '15042016';

stocks = hist_stock_data(InitialDate,EndDate,tickerMM,tickerBM,tickerEM,...
	'frequency',freq);




Returns = cell([1 length(stocks)]);
for i = 1 : length(stocks)
	Returns{i} = (stocks(i).AdjClose(2:end) - stocks(i).AdjClose(1:end-1))...
		./ stocks(i).AdjClose(1:end-1);
end
minLength = min(cellfun('length',Returns)); % get the minimum length of each cell
Returns = cellfun(@(x) x(1:minLength),Returns,'UniformOutput',false);
Returns = cell2mat(Returns);
meanReturn = mean(Returns)
yearlyReturn = (1+meanReturn).^4 - 1;
cov(Returns)