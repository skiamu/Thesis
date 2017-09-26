% run DPalgorithm (test)
clc; close all; clear variables;
%% calibration
mu = [3.38; 5; 13.03] / 100; % expected annual return
sigma = [0.31; 5.04; 10.41] / 100; % annual volatility
R = [1 0.16 -0.33; % correlation matrix
	0.16 1 0.22;
	-0.33 0.22 1];
D = diag(sigma);
S = D*R*D; % covariance matrix
param.mu = mu; param.S = S;

%% DP algo
theta = 0.05;
eta = 1e-3;
X = {1, ((1+theta):eta:1.3)',((1+theta)^2:eta:1.3)',((1+theta)^3:eta:1.3)'};
N = 3; M = 3;
tic
[ U, J] = DPalgorithm(N,M,X,param);
toc
%% plot results
figure
plot(X{2},U{2}(:,1),'r',X{2},U{2}(:,2),'b',X{2},U{2}(:,3),'y')