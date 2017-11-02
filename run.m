% run script
clc; close all; clear variables;
addpath(genpath(pwd))
rng default;
%% 1) get data from Yahoo
freqIn = 'wk'; % time-series frequency 
% freqOut = 12; % return frequency (quarterly) used in the model
% Returns = getTimeSeries(freqIn,freqOut); % get time-series from Yahoo
Returns = getTimeSeries(freqIn); % get time-series from Yahoo

%% 2) model calibration
model = 'Mixture'; % select from 'Mixture','Gaussian','GH'
GHmodel = 'NIG';
[ param, Returns, Stat, CalibrationData] = CalibrationReturns( Returns,model,GHmodel);
% save('Yahoo.mat','Returns');
debug = false;
if debug
	param = cell([2 1]);
	mu1 = [0.611;1.373;2.340]/1e3;
	mu2 = [0.683;-16.109;-17.507]/1e3;
	sigma1 = [0.069;5.666;19.121]/1e3;
	sigma2 = [0.062;6.168;52.513]/1e3;
	rho = [1 0.0633 0.0207;
		0.0633 1 -0.0236;
		0.0207 -0.0236 1];
	D1 = diag(sigma1);
	D2 = diag(sigma2);
	S1 = D1 * rho * D1;
	S2 = D2 * rho * D2;
	param{1} = {mu1,S1,0.98};
   param{2} = {mu2,S2,0.02};
end
%% 3) Dynamic Programming Algorithm
% 3.1) set parameters
N = 52; % number of time step
M = 3; % asset allocation dimension
theta = 0.07; % yearly target return
% theta = (1+theta)^(1/12)-1;
eta = 1e-4; % target set discretization
[ X ] = makeTargetSet(N,theta,eta);
VaR = 0.07; % monthly
VaR = VaR / 2 ; % weekly
alpha = 0.01;

% 3.2) run the algorithm
tic
[ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha);
toc
p_star = J{1}; % reachability problem probability

% 3.3) plot results
% I_E = 0; I_B = 0; I_C = 0;
% T = 52:(-7):1;
% N = length(T);
% s = 0;
% for yy = T
% 	s = s + 1;
% 	subplot(4,2,s);
% 	idx = find(X{yy} <= 1.14 & X{yy} >= 0.98);
% 	area(X{yy}(idx),U{yy}(idx,:))
% % 	saveas(gcf,[pwd strcat('/Latex/k',num2str(yy-1),'.png')]);
% 	title(strcat('k = ',num2str(yy-1)))
% 	% compute the areas
% 	x = X{yy}(idx);
% 	U_E = U{yy}(idx,3);
% 	U_B = U{yy}(idx,2);
% 	U_C = U{yy}(idx,1);
% 	I_E = I_E + trapz(x,U_E);
% 	I_B = I_B + trapz(x,U_B);
% 	I_C = I_C + trapz(x,U_C);
% end
save(strcat('U_',model,'.mat'),'U','p_star');
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


