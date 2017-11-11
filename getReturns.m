function [Returns,SampleStats] = getReturns( freq, M )
%getReturns is a function for downloading timeseries from Yahoo Finance and
%compute the returns for each timeseries.
%   INPUT:
%      freq = desired return frequency, freq must be one of the following
%             string {'d','wk','m','q','s','y'} [string]
%      M = number of timeseries [integer]
%   OUTPUT:
%      Returns = matrix of asset class returns in the desired frequency
%      SampleStats = 
%   REMARKS:
%      1) we initially use cellarray to store returns since timeseries may
%         be of different length

%% 1) set stocks parameters
% for more information on the timeseries check https://finance.yahoo.com
% tickerMM = 'BSV';
% tickerBM = 'BND';
tickerMM = 'SHV'; % money market
% tickerMM = 'BIL';
tickerBM = 'BGBRX'; % bond market
% tickerBM = 'FBIDX';
tickerEM = '^GSPC'; % equity market
InitialDate = '23012008';
EndDate = '15042016';
Returns = cell([1 M]); % initialization
% check input correctness
assert(any(strcmp(freq,{'d','wk','m','q','s','y'})),'Invalid input frequency : %s ',freq);

%% 2) compute Returns
if any(strcmp(freq,{'d','wk'})) % daily or weekly freq
	stocks = hist_stock_data(InitialDate,EndDate,tickerMM,tickerBM,tickerEM,...
		'frequency',freq); % data are already in the desired frequency
	for i = 1 : M
		Returns{i} = stocks(i).AdjClose(2:end) ./ stocks(i).AdjClose(1:end-1) - 1;
	end	
else % other frequencies
	freqDownload = 'wk';
	stocks = hist_stock_data(InitialDate,EndDate,tickerMM,tickerBM,tickerEM,...
		'frequency',freqDownload);
	switch freq
		case 'm' % monthly
			freqOut = 4;
		case 'q' % quarterly
			freqOut = 4 * 3;
		case 's' % semiannualy
			freqOut = 4 * 6;
		case 'y' % yearly
			freqOut = 52;
	end
	for i = 1 : M
		assert(length(stocks(i).AdjClose) >= freqOut,'There are not enough data from the freq selected');
		Returns{i} = stocks(i).AdjClose(freqOut:freqOut:end) ./ stocks(i).AdjClose(1:freqOut:end-freqOut) - 1;
	end
end

%% 3) from cell to matrix
% convert cell array to matrix: need to pay attention since the time
% series may not have the same length
minLength = min(cellfun('length',Returns)); % get the minimum length of each cell
maxLength = max(cellfun('length',Returns)); % get the maximum length of each cell
% returns are stored from left to right and we need the last minLength
% returns
Returns = cellfun(@(x) x((length(x)-minLength+1):end),Returns,'UniformOutput',false);
Returns = cell2mat(Returns); % convert cell to matrix
disp(['maximum return length = ', num2str(maxLength)]);
disp(['return length =  ',num2str(minLength)]);

%% basic Returns Stats
SampleMoments = zeros([4 M]); % each column will contain mean, std, skwe and kurtosis
SampleMoments(1,:) = mean(Returns);
SampleMoments(2,:) = std(Returns);
SampleMoments(3,:) = skewness(Returns);
SampleMoments(4,:) = kurtosis(Returns);
SampleCorr = corr(Returns);
SampleStats.Moments = SampleMoments;
SampleStats.Corr = SampleCorr;

HZmvntest(Returns); % multivariate normality test

% return distribution histogram
nbins = 15;
hist(Returns,nbins); % check return's distribution
legend('Money','Bond','Equity')
title('asset class histogram');

save('/home/andrea/Thesis/testScript/Yahoo.mat','Returns')
end % getReturns

