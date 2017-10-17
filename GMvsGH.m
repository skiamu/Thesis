% GH vs GM
clc; close all; clear variables;
addpath(genpath(pwd))

%% download data
freqIn = 'wk'; % time-series frequency 
% freqOut = 12; % return frequency (quarterly) used in the model
% Returns = getTimeSeries(freqIn,freqOut); % get time-series from Yahoo
Returns = getTimeSeries(freqIn); % get time-series from Yahoo

%% calibrate the models

model = 'Gaussian';
[ paramG, ~, StatG, CalibrationDataG] = CalibrationReturns( Returns,model);

model = 'GH'; % select from 'Mixture','Gaussian','GH'
GHmodel = 'NIG';
[ paramGH, ~, StatGH, CalibrationDataGH] = CalibrationReturns( Returns,model,GHmodel);

model = 'Mixture'; % select from 'Mixture','Gaussian','GH'
[ paramGM, Returns, StatGM, CalibrationDataGM] = CalibrationReturns( Returns,model,GHmodel);

u =[1 ; 1; 1]/3;
%% G
muG = u' * paramG.mu;
sigmaG = sqrt(u' * paramG.S * u);

%% GM
% mu must be 2x1, sigma 1x2
mu = [u' * paramGM{1}{1}; u' * paramGM{2}{1}];
sigma = zeros(1,1,2);
sigma(1,1,1) = u' * paramGM{1}{2} * u;
sigma (1,1,2) = u' * paramGM{2}{2} * u ;
p = [paramGM{1}{3} paramGM{2}{3}];

obj = gmdistribution(mu,sigma,p);
xx = (-0.1:0.0001:0.1)';
figure
plot(xx,pdf(obj,xx),xx,normpdf(xx,muG,sigmaG))
legend('Mixture','Gaussian')

%% GH

figure
plot(xx,GHdensity(paramGH,u,xx),xx,normpdf(xx,muG,sigmaG))
legend('GH','Gaussian')


figure
plot(xx,GHdensity(paramGH,u,xx),xx,normpdf(xx,muG,sigmaG),xx,pdf(obj,xx),'linewidth',1.3)
legend('GH','Gaussian','Mixture')

%% skweness kurtosis

% changing parametrization
lambda = paramGH.lambda; Chi = paramGH.Chi; Psi = paramGH.Psi;
mu = paramGH.mu; sigma = paramGH.sigma; gamma = paramGH.gamma;
% 2) compute parameters f(x,u,w(k+1)) distribution (univariate)
mu_bar = u' * mu;
sigma_bar = u' * sigma * u;
gamma_bar =  u' * gamma;

beta = gamma_bar / sigma_bar;
delta = sqrt(Chi*sigma_bar);
alpha = sqrt((Psi+gamma_bar^2 / sigma_bar)/sigma_bar);

[m, v, s, k] = nigstats(alpha, beta, mu_bar, delta)

X1 = zeros([2 1]); X2 = zeros([2 1]); X3 = zeros([2 1]);
for i = 1 : 2
	X1(i) = u' * paramGM{i}{1};
	X2(i) = u' * paramGM{i}{2} * u;
	X3(i) = paramGM{i}{3};
end
[ muX,sigma2X,gammaX,kappaX ] = ComputeMomentsGM(X1,X2,X3)













