function [param,CalibrationData] = GMcalibrationEM( Returns,k )
%GMcalibrationEM calibrates a gaussian mixture model by using the
%Expectation-Maximization algorithm
%   INPUT:
%      Returns = asset class returns [matrix]
%      k = number of mixture components
%   OUTPUT:
%      param = struct array of parameters, each component of the array is a
%              struct with the following fields: mu, S, lambda
%      CalibrationData = struct with calibration information
%   REMARKS:

Nreplicates = 150;
GM = fitgmdist(Returns,k,'Replicates',Nreplicates);
param(k) = struct(); % initialization
for i = 1 : k
	param(i).mu = GM.mu(i,:)';
	param(i).S = GM.Sigma(:,:,i);
	param(i).lambda = GM.ComponentProportion(i);
end
CalibrationData.LogL = -GM.NegativeLogLikelihood;
CalibrationData.AIC = GM.AIC;
end

