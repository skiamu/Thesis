function [param] = Calibration(MarketData,model,J_jump,dt)
%Calibration is a function for calibrating different Discrete-Event models
%   INPUT:
%      MarketData = either asset prices or asset returns [vector]
%      model = {basic,ext1,ext2} [string]
%   OUTPUT:
%      param = struct of model parameter

switch model
	case 'basic' % need to calibrate p and lambda, Market data must be prices
		[param] = CalibrationBasic(MarketData,J_jump,dt);
	case 'ext1' % need to calibrate mu and sigma, MarketData must be log-returns
		[param] = CalibrationGBM(log(1+MarketData),dt);
	case 'ext2'
		CalibrationMethod = 'LS';
		[param] = CalibrationVasicek(MarketData.LIBOR,dt,CalibrationMethod);
		[paramBasic] = CalibrationBasic(MarketData.RiskyAsset,J_jump,dt);
		param.p = paramBasic.p; param.lambda = paramBasic.lambda;
	otherwise
		error('invalid model %s',model)
		
end
end % Calibration

