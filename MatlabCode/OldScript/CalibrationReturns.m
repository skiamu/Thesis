function [ param, Returns, Stat, CalibrationData] = CalibrationReturns( Returns,model,GHmodel)
%CalibrationReturns is a function for calibrating an asset return model
%on market data
%   INPUT :
%      Returns = cell array of asset returns (X{1} = money market return, X{2} bond
%          market return, X{3} = equity market return)
%      model =  string to specify the model for the asset return:
%               'Gaussian','Mixture','GH'
%   OUTPUT :
%      param = model parameters in the format suitable for the following
%              analysis (according to the model)
%      Returns = asset class returns matrix
%      Stat = struct of time-series statistics
%      GMM = gaussian mixture object
% REMARKS:
%    1) the Return's frequency in input is the freqency specified when
%    downloading the data. Here we need to calibrate a return model with
%    frequency specified by the rebalancing policy. The two frequency may
%    disagree --> By using the compunded rule we need to convert the first
%    frequency into the second


% convert cell array to matrix: need to pay attention since the time
% series may not have the same length
minLength = min(cellfun('length',Returns)); % get the minimum length of each cell
Returns = cellfun(@(x) x(1:minLength),Returns,'UniformOutput',false);
Returns = cell2mat(Returns); % convert cell to matrix
%% model calibration
[~, d] = size(Returns);
CalibrationType = 'MM';
switch model
	case 'Mixture'
		if strcmp(CalibrationType,'MM')
% 			Mean = mean(Returns);
% 			Variance = var(Returns);
% 			Skewness = skewness(Returns);
%          Kurtosis = kurtosis(Returns);
% 			Rho = corr(Returns);
% 			Sample = zeros([15 1]);
% 			k = 0;
% 			for i = 1 : 3
% 				Sample(4*k + 1) = Mean(i);
% 				Sample(4*k + 2) =  Variance(i);
% 				Sample(4*k + 3) = Skewness(i);
% 				Sample(4*k + 4) = Kurtosis(i);
% 				k = k + 1;
% 			end
% 			Sample(13) = Rho(1,2); Sample(14) = Rho(1,3); Sample(15) = Rho(2,3);
			Sample = [0.0324;0;0;3;0.0546;0.0445;-0.46;4.25;0.1062;0.1477;-0.34;5.51;0;0;0.0342];
         [ ~, error, ~, param ] = MethodofMomentsGM( Sample,2,3,'y' );
			CalibrationData = min(error);
		else
			% 1) run calibration
			k = 2; % number of gaussian mixture components
			% 		GMM = fitgmdist(Returns,k,'Start',S); % ER algorithm
			% 		GMM = fitgmdist(Returns,k); % ER algorithm
			GMM = fitgmdist(Returns,k,'Start','plus','Replicates',20'); % ER algorithm
			% 2) build the param cell array
			param = {{GMM.mu(1,:)';GMM.Sigma(:,:,1);GMM.ComponentProportion(1)};
				{GMM.mu(2,:)';GMM.Sigma(:,:,2);GMM.ComponentProportion(2)}};
			% 3) return calibration information
			CalibrationData = GMM;
		end
	case 'GH'
		% 1) run calibration
		toll = 1e-8; maxiter = 2000; % optimization parameters
		[theta,LogL,exitFlag,numIter] = MCECMalgorithm(toll,maxiter,Returns,GHmodel);
		% 2) build the param struct
		param.lambda = theta{numIter}{1};
		param.alpha = theta{numIter}{2};
		param.mu = theta{numIter}{3};
		param.sigma = theta{numIter}{4};
		param.gamma = theta{numIter}{5};
		param.Chi = theta{numIter}{6};
		param.Psi = theta{numIter}{7};
		nParam =  1 + d + d*(d+1) / 2 + d; %NIG
		AIC = -2 * LogL + 2 * nParam;
		% 3) return calibration information
		CalibrationData = {theta,LogL,exitFlag,numIter,AIC};
	case 'Gaussian'
		param.mu = mean(Returns)';
		param.S = cov(Returns);
		CalibrationData = ecmnobj(Returns,param.mu,param.S);
	otherwise
		error('invalid model %s',model )	
end
%% some Stats
Stat.mean = mean(Returns);
Stat.Sigma = cov(Returns);
Stat.Skewness = skewness(Returns);
Stat.Kurtosis = kurtosis(Returns);

HZmvntest(Returns); % multivariate normality test

% return distribution histogram
nbins = 10;
hist(Returns,nbins); % check return's distribution
legend('Money','Bond','Equity')
title('asset class histogram');


end % function

