function [param,CalibrationData ] = modelCalibration( Returns,model,M )
%modelCalibration
%   INPUT:
%      Returns = asset class returns [matrix]
%      model = return model [string]
%   OUTPUT:
%      param = 
%      CalibrationData =

switch model
	case 'Gaussian' % gaussian model
		[param,CalibrationData] = Gcalibration(Returns);
	case 'Mixture' % gaussian mixture model
		k = 2; % number of gaussian components
		CalibrationType = 'MM';
		if strcmp(CalibrationType,'MM') % method of moments
			[param,CalibrationData] = GMcalibrationMM(Returns,k,M);
		elseif strcmp(CalibrationType,'EM') % Expectation-Maximization
			[param,CalibrationData] = GMcalibrationEM(Returns,k);
		elseif strcmp(CalibrationType,'ML') % maximum Likelihood
			[param,CalibrationData] = GMcalibrationML(Returns,k,M);
		else
			error('Invalid CalibrationType : %s', CalibrationType);
		end
	case {'GH','H','NIG','t'} % generalized Hyperbolic model
		toll = 1e-8; maxiter = 2000; % optimization parameters
		 if any(strcmp(model,{'NIG','GH','H'}))
			 [param,CalibrationData] = MCECMalgorithm(toll,maxiter,Returns,model);
		 elseif strcmp(model,'t')
			 [param, CalibrationData] = MCECMalgorithm_t(toll,maxiter,Returns,model);
		 end
	otherwise
		error('invalid model %s', model)
end

end % modelCalibration

