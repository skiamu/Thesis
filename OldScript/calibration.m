% calibration
clc; close all; clear variables;
%% market data
% import the time series and compute mean, std and correlation matrix. If
% data are not available simulate them
% ER and SD are annualized, if we want to have a quarterly rebalancing we
% need to model the distribution of quarterly return. Using the scaling
% rule we get the quarterly figures.
mu = [3.24; 5.46; 10.62]/100/4; % expected yearly return (money,bond,equity)
vol = [1e-4; 4.45; 14.77]/100/2; % volatility
R = [1 0 0; 0 1 0.0342; 0 0.0342 1]; % correlation matrix
D = diag(vol);
S = D' * R * D; % correlation matrix
N = 52 * 19; % number of weekly observations
XX = mvnrnd(mu,S,N); % data simulation
%% calibration
GMM = fitgmdist(XX,2)
