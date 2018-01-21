% GH vs GM
clc; close all; clear variables;
addpath(genpath(pwd))
rng default
%% 1) download data
M = 3;
freq = 'wk'; % return frequency, select from {'d','wk','m','q','s','y'}
[Returns,~] = getReturns( freq, M ); % download and compute asset returns
%% 2) calibrate the models

model = 'Gaussian';
[ paramG, CalibrationDataG] = modelCalibration(Returns,model,M);

model = 'NIG'; % select from 'Mixture','Gaussian','GH'
[ paramGH, CalibrationDataGH] = modelCalibration(Returns,model,M);

model = 'Mixture'; % select from 'Mixture','Gaussian','GH'
CalibrationType = 'EM';
[ paramGM, CalibrationDataGM] = modelCalibration(Returns,model,M,CalibrationType);

%% 3) Compute portfolio return
u =[0.27 ; 0.; 0.73]; % set portfolio weigths
u =[1 ; 1; 1]/3; % set portfolio weigths

x = 1;
% Gaussian model
muG = x * (1 + u' * paramG.mu);
sigmaG = sqrt(x^2 * u' * paramG.S * u);

% Mixture model
muGM = [x * (1 + u' * paramGM(1).mu); x * (1 + u' * paramGM(2).mu)]; % mu must be 2x1, sigma 1x2
sigmaGM = sqrt([x^2 * u' * paramGM(1).S * u; x^2 * u' * paramGM(2).S * u]);
p = [paramGM(1).lambda; paramGM(2).lambda];

% Generalised Hyperbolic
lambda = paramGH.lambda; Chi = paramGH.Chi; Psi = paramGH.Psi;
mu = paramGH.mu; sigma = paramGH.sigma; gamma = paramGH.gamma;
mu_bar = x * (1 + u' * mu);
sigma_bar = x^2 * u' * sigma * u;
gamma_bar =  x * u' * gamma;
% change parametrization
beta = gamma_bar / sigma_bar;
delta = sqrt(Chi*sigma_bar);
alpha = sqrt((Psi+gamma_bar^2 / sigma_bar)/sigma_bar);

% estimated moments
PortfolioReturns = x * (1 + Returns * u);
EstimatedSkw = skewness(PortfolioReturns);
EstimatedKurt = kurtosis(PortfolioReturns);
EstimatedMean = mean(PortfolioReturns);
EstimatedVar = var(PortfolioReturns);
%% 4) Quantitative analysis
% moments NIG
[meanNIG, varNIG, skwNIG, kurtNIG] = nigstats(alpha, beta, mu_bar, delta);
% moments GM
[ meanGM,volGM,skwGM,kurtGM] = ComputeMomentsGM(muGM,sigmaGM,p);

disp('%%%%%%%%%%%%%%% NIG %%%%%%%%%%%%%%%%%%%%%%')
disp(['meanNIG = ',num2str(meanNIG)]);
disp(['volNIG = ',num2str(sqrt(varNIG))]);
disp(['skwenessNIG = ',num2str(skwNIG)]);
disp(['kurtosisNIG = ',num2str(kurtNIG)]);
disp(['logL = ',num2str(CalibrationDataGH.LogL)]);
disp(['AIC = ',num2str(CalibrationDataGH.AIC)]);
disp('%%%%%%%%%%%% Gaussian MIxture %%%%%%%%%%%%')
disp(['meanGM = ',num2str(meanGM)]);
disp(['volGM = ',num2str(volGM)]);
disp(['skwenessGM = ',num2str(skwGM)]);
disp(['kurtosisGM = ',num2str(kurtGM)]);
disp(['logL = ',num2str(CalibrationDataGM.LogL)]);
disp(['AIC = ',num2str(CalibrationDataGM.AIC)]);
disp('%%%%%%%%%%%%%% Gaussian %%%%%%%%%%%%%%%%%%%')
disp(['logL = ',num2str(CalibrationDataG.LogL)]);
disp('%%%%%%%%%%%% Estimated Moments %%%%%%%%%%%%%')
disp(['Estimated Skweness = ',num2str(EstimatedSkw)]);
disp(['Estimated Kurtosis = ',num2str(EstimatedKurt)]);
disp(['Estimated Vol = ',num2str(sqrt(EstimatedVar))]);
disp(['Estimated Mean = ',num2str(EstimatedMean)]);
