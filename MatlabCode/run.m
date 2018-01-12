% run file
clc; close all; clear variables;
addpath(genpath(pwd))
DefaultPlotting
% rng default;
M = 3; % asset allocation dimension
%% 1) get data from Yahoo
freq = 'wk'; % return frequency, select from {'d','wk','m','q','s','y'}
[Returns,SampleStats] = getReturns( freq, M ); % download and compute asset returns

%% 2) model Calibration
model = 'Gaussian'; % select from {'Gaussian','Mixture','GH'}
CalibrationType = 'EM'; % select from {'MM','ML','EM'} (only for GM model)
[param,CalibrationData] = modelCalibration( Returns,model,M,CalibrationType);

%% 3) Dynamic Programming Algorithm
% 3.1) set parameters
VaR = 0.07; % monthly
alpha = 0.01; % confidence level VaR
switch freq % number of time step for a 2-year investment
	case 'wk'
		N = 104;
		NstepPlot = 26;
		VaR = VaR / 2;
	case 'm'
		N = 24;
		NstepPlot = 6;
	case 'q'
		N = 12;
		NstepPlot = 3;
		VaR = VaR * sqrt(3);
end
theta = 0.07; % yearly target return
eta = 1e-3/5; % target set discretization
[ X ] = makeTargetSet(N,theta,eta);
% 3.2) run the algorithm
tic
[ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha);
toc
p_star = J{1}; % reachability problem probability

%% 4) Validation
Nsim = 1e+6;
[ w ] = SimulationReturns(param,Nsim,M,N,model);
[ p_starMC ] = Validation(w,X,U,Nsim,M,N);
SaveFlag = true;
if SaveFlag
	save(strcat('/home/andrea/Thesis/MatlabCode/',model,freq,'.mat'),'J', 'U', 'X',...
		'p_star','p_starMC')
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


