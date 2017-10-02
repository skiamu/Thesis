% run script
clc; close all; clear variables;

%% get data from Yahoo
freqIn = 'wk'; % time-series frequency 
freqOut = 12; % return frequency (quarterly) used in the model
Returns = getTimeSeries(freqIn,freqOut); % get time-series from Yahoo

%% model calibration
model = 'Mixture';
CalibrationType = 'EM';
[ param, Returns, Stat] = CalibrationReturns( Returns, CalibrationType,...
	model);

%% Dynamic Programming Algorithm
N = 8; % number of time step
M = 3; % asset allocation dimension
theta = 0.07; % yearly target return
eta = 1e-3; % target set discretization
X = cell([N+1 1]); % initialization
for i = 2 : N
	X{i} = (0.5:eta:1.9)';
end
X{1} = 1; X{end} = ((1+theta)^2:eta:1.9)';
VaR = 0.07; % monthly
VaR = VaR * sqrt(3); % quarterly
alpha = 0.01;
tic
[ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha);
toc
J{1} % reachability problem probability
%% plot results

for yy = N :-1: 2
	figure
	idx = find(X{yy} <= 1.2 & X{yy} >= 0.95);
	area(X{yy}(idx),U{yy}(idx,:))
	saveas(gcf,[pwd strcat('/Latex/k',num2str(yy-1),'.png')]);
	title(strcat('k = ',num2str(yy-1)))
end










