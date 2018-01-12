% model calibration comparison.
% This scripts calibrate a GM model to data with a known GM distribution to
% see which calibration method performs better.

clc; close all; clear variables;
addpath(genpath(pwd))
rng('default')
% rng default;
M = 3; % asset allocation dimension

%% simulate GM data
% the data is taken from Pola's paper
mu1 = [.000611 .001373 .002340];
mu2 = [.000683 -.016109 -.017507];
sigma1 = [.000069 .005666 .019121]';
sigma2 = [.000062 .006168 .052513]';
R = [1 .0633 .0207; .0633 1 -.0236; .0207 -.0236 1];
Sigma1 = corr2cov(sigma1,R);
Sigma2 = corr2cov(sigma2,R);
lambda = 0.98;
X1 = mvnrnd(mu1,Sigma1,ceil(lambda * 1e+4));
X2 = mvnrnd(mu2,Sigma2,ceil((1-lambda) * 1e+4));
X = [X1; X2];
%% model calibration
model = 'Mixture'; % select from {'Gaussian','Mixture','GH'}
CalibrationType = 'EM'; % select from {'MM','ML','EM'} (only for GM model)
[param,CalibrationData] = modelCalibration( X,model,M,CalibrationType);

ErrorMu1 = abs((param(1).mu - mu1') ./ mu1' * 100);
ErrorMu2 = abs((param(2).mu - mu2') ./ mu2' * 100);
ErrorSigma1 = abs((param(1).S - Sigma1) ./ Sigma1 * 100);
ErrorSigma2 = abs((param(2).S - Sigma2) ./ Sigma2 * 100);
ErrorLambda = abs(param(1).lambda - lambda) / lambda * 100;




