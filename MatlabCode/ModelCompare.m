% GH vs GM
clc; close all; clear variables;
addpath(genpath(pwd))
% rng default
%% 1) download data
freq = 'wk'; % return frequency, select from {'d','wk','m','q','s','y'}
[Returns,~] = getReturns( freq, M ); % download and compute asset returns

%% 2) calibrate the models

model = 'Gaussian';
[ paramG, ~, StatG, CalibrationDataG] = CalibrationReturns( Returns,model);

model = 'GH'; % select from 'Mixture','Gaussian','GH'
GHmodel = 'NIG';
[ paramGH, ~, StatGH, CalibrationDataGH] = CalibrationReturns( Returns,model,GHmodel);

model = 'Mixture'; % select from 'Mixture','Gaussian','GH'
[ paramGM, Returns, StatGM, CalibrationDataGM] = CalibrationReturns( Returns,model,GHmodel);

save('Yahoo.mat','Returns');
nbins = 10;
hist(Returns,nbins); % check return's distribution
legend('Money','Bond','Equity')
title(strcat('asset class histogram,','freq = ',freqIn));
% saveas(gcf,[pwd strcat('/Latex/secondWIP/histReturns',freqIn,'.png')]);
%% 3) Compute portfolio return
u =[0.0 ; 0.27; 0.73]; % set portfolio weigths
x = 1.097;
% Gaussian model
muG = x * (1 + u' * paramG.mu);
sigmaG = sqrt(x^2 * u' * paramG.S * u);

% Mixture model
muGM = [x * (1 + u' * paramGM{1}{1}); x * (1 + u' * paramGM{2}{1})]; % mu must be 2x1, sigma 1x2
sigmaGM = zeros(1,1,2);
sigmaGM(1,1,1) = x^2 * u' * paramGM{1}{2} * u;
sigmaGM (1,1,2) = x^2 * u' * paramGM{2}{2} * u ;
p = [paramGM{1}{3} paramGM{2}{3}];

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
%% 4) Graphical analysis
% 4.1) G vs GM
GMobj = gmdistribution(muGM,sigmaGM,p); % create a GaussianMixture object
xx = (1.05:0.0001:1.15)';
figure
plot(xx,pdf(GMobj,xx),'.',xx,normpdf(xx,muG,sigmaG),'linewidth',1.3)
h = legend('Mixture','Gaussian','Location','best'); h.FontSize = 14;
title(['G vs GM, freq = ',freqIn])
xlabel('portfolio value')
grid on
% 4.2) G vs GH
figure
plot(xx,GHdensityPortfolioReturns(paramGH,u,xx,x),'.',xx,normpdf(xx,muG,sigmaG),...
	'linewidth',1.3)
h = legend('GH','Gaussian'); h.FontSize = 14;
title(['G vs GH, freq = ',freqIn])
xlabel('portfolio value')
grid on
% 4.3) G vs GH vs GM
figure
plot(xx,GHdensityPortfolioReturns(paramGH,u,xx,x),'.',xx,normpdf(xx,muG,sigmaG),...
	xx,pdf(GMobj,xx),'.','linewidth',1.3)
h = legend('GH','Gaussian','Mixture','Location','NorthEast'); h.FontSize = 14;
title(['G vs GM vs GH, freq = ',freqIn])
xlabel('portfolio value')
grid on
% saveas(gcf,[pwd strcat('/Latex/secondWIP/GvsGMvsGH',freqIn,'.png')]);
% 4.4) Histogram
figure
histogram(PortfolioReturns,140)
hold on
plot(xx,GHdensityPortfolioReturns(paramGH,u,xx,x),'.',xx,normpdf(xx,muG,sigmaG),...
	xx,pdf(GMobj,xx),'.','linewidth',1.3)
h = legend('Hist','GH','Gaussian','Mixture'); h.FontSize = 14;
title('G vs GH vs GM')
xlabel('portfolio value')
grid on
% saveas(gcf,[pwd strcat('/Latex/secondWIP/hist',freqIn,'.png')]);

% 4.5) 
% u1 = [0.27; 0; 0.73]; % cash
% u2 = [0; 0.27; 0.73]; % bond
% muGM = [x * (1 + u1' * paramGM{1}{1}); x * (1 + u1' * paramGM{2}{1})]; % mu must be 2x1, sigma 1x2
% sigmaGM = zeros(1,1,2);
% sigmaGM(1,1,1) = x^2 * u1' * paramGM{1}{2} * u1;
% sigmaGM (1,1,2) = x^2 * u1' * paramGM{2}{2} * u1 ;
% figure
% plot(xx,pdf(GMobj,xx),'-','linewidth',1.3);
% hold on
% muGM = [x * (1 + u2' * paramGM{1}{1}); x * (1 + u2' * paramGM{2}{1})]; % mu must be 2x1, sigma 1x2
% sigmaGM = zeros(1,1,2);
% sigmaGM(1,1,1) = x^2 * u2' * paramGM{1}{2} * u2;
% sigmaGM (1,1,2) = x^2 * u2' * paramGM{2}{2} * u2 ;
% plot(xx,pdf(GMobj,xx),'-','linewidth',1.3);
% legend('Cash','Bond')
%% 5) Quantitative analysis
% moments NIG
[meanNIG, varNIG, skwNIG, kurtNIG] = nigstats(alpha, beta, mu_bar, delta);
% moments GM
[ meanMG,varGM,skwGM,kurtGM] = ComputeMomentsGM(paramGM,u,x);

disp('%%%%%%%%%%%%%%% NIG %%%%%%%%%%%%%%%%%%%%%%')
disp(['skwenessNIG = ',num2str(skwNIG)]);
disp(['kurtosisNIG = ',num2str(kurtNIG)]);
disp(['logL = ',num2str(CalibrationDataGH{2})]);
disp(['AIC = ',num2str(CalibrationDataGH{5})]);
disp('%%%%%%%%%%%% Gaussian MIxture %%%%%%%%%%%%')
disp(['skwenessGM = ',num2str(skwGM)]);
disp(['kurtosisGM = ',num2str(kurtGM)]);
disp(['logL = ',num2str(CalibrationDataGM.NegativeLogLikelihood)]);
disp(['AIC = ',num2str(CalibrationDataGM.AIC)]);
disp('%%%%%%%%%%%%%% Gaussian %%%%%%%%%%%%%%%%%%%')
disp(['logL = ',num2str(CalibrationDataG(1))]);
disp('%%%%%%%%%%%% Estimated Moments %%%%%%%%%%%%%')
disp(['Estimated Skweness = ',num2str(EstimatedSkw)]);
disp(['Estimated Kurtosis = ',num2str(EstimatedKurt)]);
disp(['Estimated Variance = ',num2str(EstimatedVar)]);
disp(['Estimated Mean = ',num2str(EstimatedMean)]);
