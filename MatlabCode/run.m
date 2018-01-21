% run file
clc; close all; clear variables;
addpath(genpath(pwd))
DefaultPlotting
rng default;
M = 3; % asset allocation dimension
%% 1) get data from Yahoo
freq = 'wk'; % return frequency, select from {'d','wk','m','q','s','y'}
[Returns,SampleStats] = getReturns( freq, M ); % download and compute asset returns

%% 2) model Calibration
model = 'Mixture'; % select from {'Gaussian','Mixture','GH'}
CalibrationType = 'EM'; % select from {'MM','ML','EM'} (only for GM model)
[param,CalibrationData] = modelCalibration( Returns,model,M,CalibrationType);

% horizon correction
freq = 'm';
[param] = HorizonCorrection(freq,param,model);
%% 3) Dynamic Programming Algorithm
% 3.1) set parameters
SetODAAparameters
% 3.2) run the algorithm
tic
[ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha);
time = toc;
p_star = J{1}; % reachability problem probability

%% 4) Validation
Nsim = 1e+5; % numero of montecarlo simulation for every time period
[ w ] = SimulationReturns(param,Nsim,M,N,model);
[ p_starMC,Statistics] = Validation(w,X,U,Nsim,M,N,freq);
SaveFlag = true;
if SaveFlag
	save(strcat('/home/andrea/Thesis/MatlabCode/',model,freq,'.mat'),'J', 'U', 'X',...
		'p_star','p_starMC','time','Statistics');
end
%% 5) plot results
k = 0;
for i = N : -NstepPlot : 1
	k = k + 1;
	subplot(3,2,7-k)
	idx = find(X{i} <= 1.3 & X{i} >= 0.90);
	area(X{i}(idx),U{i}(idx,:))
	grid on
	ylim([0 1])
	xlim([0.9 1.3])
	xlabel('portfolio value')
	title(strcat('k = ',num2str(i-1)))
	% 	saveas(gcf,[pwd strcat('/Latex/thirdWIP/k',num2str(i-1),model,'.png')]);
end
k = k + 1;
subplot(3,2,7-k)
b = bar(diag(U{1}));
grid on
legend(b,'Cash','Bond','Equity','Location','NorthWest')
print(strcat('/home/andrea/Thesis/Latex/final/Images/maps',model,freq),...
	'-dpng', '-r900');


