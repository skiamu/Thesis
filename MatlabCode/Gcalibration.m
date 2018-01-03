function [param,CalibrationData] = Gcalibration(Returns)
%Gcalibration calibrates the gaussian model
%   INPUT:
%      Returns = asset class returns [matrix]
%   OUTPUT:
%      param = struct array of parameters, each component of the array is a
%              struct with the following fields: mu, S, lambda
%      CalibrationData = struct with calibration information

param.mu = mean(Returns)';
param.S = cov(Returns);
CalibrationData.LogL = -ecmnobj(Returns,param.mu,param.S);
end

