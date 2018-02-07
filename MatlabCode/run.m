% run file
clc; close all; clear variables;
addpath(genpath(pwd))
DefaultPlotting
rng default;
M = 3; % asset allocation dimension
%% 1) get data from Yahoo
freq = 'wk'; % return frequency, select from {'d','wk','m','q','s','y'}
[Returns,SampleStats] = getReturns( freq, M ); % download and compute asset returns
save('/home/andrea/Thesis/PythonCode/AssetClassReturn.mat','Returns')
%% 2) model Calibration
model = 'Mixture'; % select from {'Gaussian','Mixture','GH'}
CalibrationType = 'EM'; % select from {'MM','ML','EM'} (only for GM model)
[param,CalibrationData] = modelCalibration( Returns,model,M,CalibrationType);

% horizon correction
freq = 'wk';
[param,t] = HorizonCorrection(freq,param,model);
%% 3) ODAA Algorithm
% 3.1) set parameters
SetODAAparameters
% 3.2) run the ODAA algorithm
tic
[ U, J] = ODAAalgorithm(N,M,X,param,model,VaR,alpha);
time = toc;
p_star = J{1}; % reachability problem probability

%% 4) CPPI and Constant-Mix
% 4.1) synthetize Constant-MIx strategy
[U_ConstantMix] = ConstantMix(param,model,VaR,M,N,alpha);
% 4.2) synthetize CPPI strategy
ER = (1 + mean(Returns)).^t - 1;
r = ER(1); % cash return in the rebalancing frequency
u0 = U{1}'; % initial asset allocation from ODAA
m = 6; % CPPI multiplier
[U_CPPI,~,~] = CPPI(u0,X,r,m,N,param,model,VaR,alpha);
save(strcat('/home/andrea/Thesis/MatlabCode/CPPI',model,freq,'.mat'), 'U_CPPI');
%% 5) Validation Monte-Carlo
rng default
Nsim = 2e+5; % numero of montecarlo simulation for every time period
% 5.1) simulate asset class returns
[w] = SimulationReturns(param,Nsim,M,N,model);
% 5.2) validation ODAA strategy
[p_starMC,~] = Validation(w,X,U,Nsim,M,N,freq,'ODAA',r);
% 5.3) validation Constant-Mix strategy
[p_starMC_ConstantMix,~] = Validation(w,X,U_ConstantMix,Nsim,M,N,freq,'ConstantMix',r);
% 5.4) validation CPPI strategy
[p_starMC_CPPI,~] = Validation(w,X,U_CPPI,Nsim,M,N,freq,'CPPI',r);
legend('ODAA','Constant-mix','CPPI');
title('Investement return empirical density function');
print(strcat('/home/andrea/Thesis/Latex/final/Images/Densities',model,freq),...
	'-dpng', '-r900');
SaveFlag = true;
if SaveFlag
	save(strcat('/home/andrea/Thesis/MatlabCode/',model,freq,'.mat'),'J', 'U', 'X',...
		'p_star','p_starMC','time','Statistics');
end

%% 6) plot maps
k = 0;
for i = N : -NstepPlot : 1
	k = k + 1;
	subplot(3,2,7-k)
	idx = find(X{i} <= 1.3 & X{i} >= 0.60);
	area(X{i}(idx),U_CPPI{i}(idx,:))
	grid on
	ylim([0 1])
	xlim([0.7 1.3])
	xlabel('portfolio value')
	title(strcat('k = ',num2str(i-1)))
	% 	saveas(gcf,[pwd strcat('/Latex/thirdWIP/k',num2str(i-1),model,'.png')]);
end
k = k + 1;
subplot(3,2,7-k)
b = bar(diag(U{1}));
grid on
legend(b,'Cash','Bond','Equity','Location','NorthWest')
print(strcat('/home/andrea/Thesis/Latex/final/Images/mapsCPPI',model,freq),...
	'-dpng', '-r900');


