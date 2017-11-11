% run file
clc; close all; clear variables;
addpath(genpath(pwd))
% rng default;
M = 3; % asset allocation dimension
%% 1) get data from Yahoo
freq = 'wk'; % return frequency, select from {'d','wk','m','q','s','y'}
[Returns,SampleStats] = getReturns( freq, M ); % download and compute asset returns

%% 2) model Calibration
model = 'Mixture'; % select from {'Gaussian','Mixture','GH'}
[param,CalibrationData] = modelCalibration( Returns,model,M );

%% 3) Dynamic Programming Algorithm
% 3.1) set parameters
N = 24; % number of time step
theta = 0.07; % yearly target return
% theta = (1+theta)^(1/12)-1;
eta = 1e-3 / 2; % target set discretization
[ X ] = makeTargetSet(N,theta,eta);
VaR = 0.07; % monthly
VaR = VaR / sqrt(4);
alpha = 0.01;
% 3.2) run the algorithm
tic
[ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha);
toc
p_star = J{1}; % reachability problem probability

%% 4) Validation
Nsim = 1e+6;
[ w ] = SimulationReturns(param,Nsim,M,N,model);
[ p_starMC ] = Validation(w,X,U,Nsim,M,N);









