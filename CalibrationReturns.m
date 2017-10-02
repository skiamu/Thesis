function [ param, Returns, Stat] = CalibrationReturns( Returns, CalibrationType,... 
	model)
%CalibrationReturns is a function for calibrating an asset return model
%on market data
%   INPUT :
%      Returns = cell array of asset returns (X{1} = money market return, X{2} bond
%          market return, X{3} = equity market return)
%      CalibrationType =  'ER ': Expectation Maximization
%                         'MM' : Method of Moment
%      model =  string to specify the model for the asset return:
%               'Gaussian','Mixture'
%   OUTPUT :
%      param = model parameters in the format suitable for the following
%              analysis (according to the model)
%      Returns = asset class returns matrix
%      Stat = struct of time-series statistics
%
% REMARKS:
%    1) the Return's frequency in input is the freqency specified when
%    downloading the data. Here we need to calibrate a return model with
%    frequency specified by the rebalancing policy. The two frequency may
%    disagree --> By using the compunded rule we need to convert the first
%    frequency into the second

% check for impossible calibrations
if strcmp(model,'Gaussian') && strcmp(CalibrationType,'EM')
	error('invalid model and CalibrationType input');
end
	
% convert cell array to matrix: need to pay attention since the time
% series may not have the same length
minLength = min(cellfun('length',Returns)); % get the minimum length of each cell
Returns = cellfun(@(x) x(1:minLength),Returns,'UniformOutput',false);
Returns = cell2mat(Returns); % convert cell to matrix
%% model calibration
switch CalibrationType
	case 'EM'
		k = 2; % number of gaussian mixture components
		GMM = fitgmdist(Returns,k); % ER algorithm
		param = {{GMM.mu(1,:)';GMM.Sigma(:,:,1);GMM.ComponentProportion(1)};
			{GMM.mu(2,:)';GMM.Sigma(:,:,2);GMM.ComponentProportion(2)}};
	case 'MM'
		switch model
			case 'Gaussian'
				param.mu = mean(Returns)';
				param.S = cov(Returns);
			case 'Mixture'
				
			otherwise
				error('invalid model %s',model);
		end
	otherwise
		error('invalid CalibrationType %s',CalibrationType);
end

%% some Stats
Stat.mean = mean(Returns);
Stat.Sigma = cov(Returns);
Stat.Skewness = skewness(Returns);
Stat.Kurtosis = kurtosis(Returns);

cd HZmvntest;
HZmvntest(Returns); % multivariate normality test
cd ..;

nbins = 10;
hist(Returns,nbins); % check return's distribution
legend('Money','Bond','Equity')
title('asset class histogram');


end % function

