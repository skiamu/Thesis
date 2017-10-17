% run script
clc; close all; clear variables;
addpath(genpath(pwd))
%% 1) get data from Yahoo
freqIn = 'wk'; % time-series frequency 
freqOut = 12; % return frequency (quarterly) used in the model
Returns = getTimeSeries(freqIn,freqOut); % get time-series from Yahoo
% Returns = getTimeSeries(freqIn); % get time-series from Yahoo

%% 2) model calibration
model = 'Mixture'; % select from 'Mixture','Gaussian','GH'
GHmodel = 'NIG';
[ param, Returns, Stat, CalibrationData] = CalibrationReturns( Returns,model,GHmodel);
% save('Yahoo.mat','Returns');
%% 3) Dynamic Programming Algorithm
% 3.1) set parameters
N = 8; % number of time step
M = 3; % asset allocation dimension
theta = 0.07; % yearly target return
eta = 1e-3; % target set discretization
[ X ] = makeTargetSet(N,theta,eta);
VaR = 0.07; % monthly
VaR = VaR*sqrt(3); % quarterly
alpha = 0.01;

% 3.2) run the algorithm
tic
[ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha);
toc
p_star = J{1}; % reachability problem probability

% 3.3) plot results
I_E = 0; I_B = 0; I_C = 0;
for yy = N :-1: 2
	subplot(4,2,yy-1);
	idx = find(X{yy} <= 1.2 & X{yy} >= 0.95);
	area(X{yy}(idx),U{yy}(idx,:))
	saveas(gcf,[pwd strcat('/Latex/k',num2str(yy-1),'.png')]);
	title(strcat('k = ',num2str(yy-1)))
	% compute the areas
	x = X{yy}(idx);
	U_E = U{yy}(idx,3);
	U_B = U{yy}(idx,2);
	U_C = U{yy}(idx,1);
	I_E = I_E + trapz(x,U_E);
	I_B = I_B + trapz(x,U_B);
	I_C = I_C + trapz(x,U_C);
end

%% 4) Validation ODA
Nsim = 1e6; 
simulationMethod = 'built-in';
GM = CalibrationData;
[ w ] = SimulationReturns(param,Nsim,M,N,model,simulationMethod,GM );
% [ w ] = SimulationReturns(param,Nsim,M,N,model);
[p_star_MC] = Validation(w,X,U,Nsim,M,N);

%% 5) Comparison CPPI
u0 = U{1}';
MeanReturns = mean(Returns);
r = MeanReturns(1); % quarterly expected cash returns
m = 6;
[Ucppi,Floor,Cushion] = CPPI(u0,X,r,m,N,param,model,VaR,alpha);
[p_star_cppi] = Validation(w,X,Ucppi,Nsim,M,N);


