function Returns = getTimeSeries(freqIn, freqOut)
%getTimeSeries is a function for getting time series from Yahoo finance.
% In the following script we retrive financial time-series from the info
% provider Yahoo Finance. We take 3 indexes that represents the following
% asset class: 1) Money Market (iShares Short Treasury Bond ETF)
%              2) Bond Market (iShares 20+ Year Treasury Bond ETF)
%              3) Equity Market (S&P 500)
% Finally we compute the return according to the specified frequency
%
% INPUT :
%    freqIn = price frequency downloaded data(must be either'd' for daily, 'wk' for weekly, or 'mo' for monthly)
%    freqOut = return frequency in output (4 = monthy, 12 = quarterly, 24 = semiannualy, 52 = yearly)
% OUTPUT:
%    Returns = cell array of returns

%% download the data
% freq = 'wk';
tickerMM = 'SHV';
tickerBM = 'BGBRX';
tickerEM = '^GSPC';
InitialDate = '23012008';
EndDate = '15042016';

cd hist_stock_data; % access the subfolder where the function is stored (Linux command)
stocks = hist_stock_data(InitialDate,EndDate,tickerMM,tickerBM,tickerEM,...
	'frequency',freqIn);
cd ..; % Linux command

%% compute the return
% Return is a cell array: Return{1} = Money market returns
%                         Return{2} = Bond market returns
%                         Return{3} = Equity market returns
% REMARK :
% 1) we use cell array because in general the timeseries have not the
% same length (e.g. quotes not available from InitialDate for a certain stock )
% 2) Returns are computed according to the specified frequency
Returns = cell([1 length(stocks)]);
if nargin < 2 % weekly returns
	for i = 1 : length(stocks)
		Returns{i} = (stocks(i).AdjClose(2:end) - stocks(i).AdjClose(1:end-1))...
			./ stocks(i).AdjClose(1:end-1);
	end
else
	% return with specified frequency (not overlapping)
	for i = 1 : length(stocks)
		Returns{i} = (stocks(i).AdjClose(freqOut:freqOut:end) - stocks(i).AdjClose(1:freqOut:end-freqOut))...
			./ stocks(i).AdjClose(1:freqOut:end-freqOut);
	end
end

end % function