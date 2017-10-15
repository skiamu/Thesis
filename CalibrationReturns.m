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

switch model
	case 'Mixture'
		% 0) run method of moments for an initial estimate
% 		[x0] = MethodMomentsGM(Returns);
% 		D1 = diag(x0(5:7)); D2 = diag(x0(11:13));
% 		rho12 = x0(14); rho13 = x0(15); rho23 = x0(16);
% 		R = eye(3);
% 		R(1,2) = rho12; R(2,1) = rho12; R(1,3) = rho13; R(3,1) = rho13;
% 		R(2,3) = rho23; R(3,2) = rho23;
% 		Sigma = zeros([3 3 2]);
% 		Sigma(:,:,1) = D1' * R * D1; Sigma(:,:,2) = D2' * R * D2;
% 		Mu = [x0(2:4)'; x0(8:10)'];
% 		PComponents = [x0(1) 1-x0(1)];
% 		S = struct('mu',Mu,'Sigma',Sigma,'ComponentProportion',PComponents);
		
		% 1) run calibration
		k = 2; % number of gaussian mixture components
% 		GMM = fitgmdist(Returns,k,'Start',S); % ER algorithm
		GMM = fitgmdist(Returns,k,'Replicates',50); % ER algorithm
		% 2) build the param cell array
		param = {{GMM.mu(1,:)';GMM.Sigma(:,:,1);GMM.ComponentProportion(1)};
			{GMM.mu(2,:)';GMM.Sigma(:,:,2);GMM.ComponentProportion(2)}};
		% 3) return calibration information
		CalibrationData = GMM;
	case 'GH'
		% 1) run calibration
		toll = 1e-3; maxiter = 500; % optimization parameters
		[theta,LogL,exitFlag,numIter] = MCECMalgorithm(toll,maxiter,Returns,GHmodel);
		% 2) build the param struct
		param.lambda = theta{numIter}{1};
		param.Chi = theta{numIter}{2};
		param.Psi = theta{numIter}{3};
		param.mu = theta{numIter}{4};
		param.sigma = theta{numIter}{5};
		param.gamma = theta{numIter}{6};
		% 3) return calibration information
		CalibrationData = {theta,LogL,exitFlag,numIter};
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

