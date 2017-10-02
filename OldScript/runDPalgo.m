% run DPalgorithm (test)
clc; close all; clear variables;
rng(1992)
%% calibration gaussian
% mu = [3.38; 5; 13.03] / 100; % expected annual return
% sigma = [0.31; 5.04; 10.41] / 100; % annual volatility
% R = [1 0.16 -0.33; % correlation matrix
% 	0.16 1 0.22;
% 	-0.33 0.22 1];
% D = diag(sigma);
% S = D*R*D; % covariance matrix
% param.mu = mu; param.S = S;

%% calibration Mixture
% sigma1 = [0.000069; 0.00566; 0.002340]*sqrt(52);
% sigma2 = [0.000062; 0.006168; 0.052513]*sqrt(52);
% D1 = diag(sigma1); D2 = diag(sigma2);
% R = [1 0.0633 0.0207; 0.0633 1 -0.0236; 0.0207 -0.0236 1];
% S1 = D1'*R*D1; S2 = D2'*R*D2;
% param = {{[0.000611;0.001373;0.002340]*52, S1, 0.98};
% 	{[0.000683;-0.016109;-0.017507]*52, S2, 0.02}};
%% calibration
% import the time series and compute mean, std and correlation matrix. If
% data are not available simulate them
% ER and SD are annualized, if we want to have a quarterly rebalancing we
% need to model the distribution of quarterly return. Using the scaling
% rule we get the quarterly figures.
% data from paper : non-gaussian world
mu = [3.24; 5.46; 10.62]/100/4; % expected yearly return (money,bond,equity)
vol = [0.5; 4.45; 14.77]/100/2; % volatility
R = [1 0 0; 0 1 0.0342; 0 0.0342 1]; % correlation matrix
D = diag(vol);
S = D' * R * D; % correlation matrix
N = 52 * 19; % number of weekly observationsGMM = fitgmdist(XX,2);

XX = mvnrnd(mu,S,N); % data simulation
XX_T = mvtrnd(R,3,N) + ones([N 1]) * mu';
GMM = fitgmdist(XX_T,2);
model = 'Mixture';

%% DP algo
N = 8; M = 3;
theta = 0.08;
eta = 1e-3;
X = cell([N+1 1]); % initialization
for i = 2 : N
	X{i} = (0.5:eta:1.9)';
end
X{1} = 1; X{end} = ((1+theta)^2:eta:1.9)';
if strcmp(model,'Mixture')
	param = {{GMM.mu(1,:)';GMM.Sigma(:,:,1);GMM.ComponentProportion(1)};
	{GMM.mu(2,:)';GMM.Sigma(:,:,2);GMM.ComponentProportion(2)}};
elseif strcmp(model,'Gaussian')
	 param.mu = mu; param.S = S;
end
% X = {1, ((1+theta):eta:1.9)',((1+theta)^2:eta:1.9)',((1+theta)^3:eta:1.9)'};
% X = {1, (0.9:eta:1.3)',(0.9:eta:1.3)',((1+theta)^3:eta:1.3)'};
tic
[ U, J] = DPalgorithm(N,M,X,param,model);
toc
%% plot results
figure
% plot(X{2},U{2}(:,1),'r',X{2},U{2}(:,2),'b',X{2},U{2}(:,3),'y')
yy = 6;
idx = find(X{yy} <= 1.2 & X{yy} >= 0.95);
area(X{yy}(idx),U{yy}(idx,:))
J{1}







