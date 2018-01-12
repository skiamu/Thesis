% test script
clc; close all; clear variables;
addpath(genpath(pwd))
load('LIBOR.mat')
load('SP500Fut.mat')
idx1 = find(DateSP500Fut == '22-Jan-2010');
idx2 = find(DateSP500Fut == '25-Apr-2016');
S = SP500Fut(idx1:idx2);
Returns = S(2:end) ./ S(1:end-1) - 1;
%% discrete dynamics
model = 'basic'; % select from {'basic','ext1','ext2'}
% freq = 'd';
% M = 3;
% [Returns,~,stocks] = getReturns( freq, M );
% S = stocks(3).AdjClose;
% save('/home/andrea/Thesis/PythonCode/RiskyAsset.mat','S');
J_jump = 0.07; % jump size
dt = 1 / 252; % time in years
switch model
	case 'basic'
		[param] = Calibration(S,model,J_jump,dt);
		p = param.p; lambda = param.lambda;param.r = 0.055;
	case 'ext1'
		[param] = Calibration(Returns,model,J_jump,dt);
		param.r = 0.055; % cash return
		mu = param.mu; sigma = param.sigma;
		mu_tilde = mu - 0.5 * sigma^2;
		param.p = (exp(2 * mu_tilde * J_jump / sigma^2) - 1 ) / ...
			(2 * sinh(2 * mu_tilde * J_jump / sigma^2)); % probability positive jump
	case 'ext2'
		MarketData.LIBOR = LIBOR(1:2000) / 100;
		MarketData.RiskyAsset = S;
		[param] = Calibration(MarketData,model,J_jump,dt);
end

%% DP algorithm
N = 10; % number of events
theta = 0.07; % yearly target return
eta = 1e-2/5; % target set discretization
[ X ] = makeTargetSet(N,theta,eta);
tic
[U,J] = DPalgorithmDES(N,X,J_jump,param,model);
toc

for k = 2 : 10
	idx = find(X{k} <= 1.8 & X{k} >= 0.5);
	subplot(5,2,k)
	area(X{k}(idx),U{k}(idx))
	title(strcat('k = ',num2str(k-1)))
	grid on
end
saveas(gcf,'/home/andrea/Thesis/Latex/thirdWIP/mapsBasic.png');
p_star = J{1};
% save variables J, U, X and p_star
save('/home/andrea/Thesis/MatlabCode/DiscreteEvent/basic.mat','J', 'U', 'X', 'p_star')
Nsim = 1e+6;
[p_starMC,Times] = ValidationDES(X,U,Nsim,N,param,J_jump);
%% plot
u = 0.8;
x = 1.2;
csi = x * J_jump * u;
z = (0.5:0.0001:1.9)';
plot(z,pfDESext1(z,x,u,J_jump,param),'r.-')
tic
trapz(z,pfDESext1(z,x,u,J_jump,param))
toc
% simpsons(@(z)pfDESext2(z,x,u,J_jump,param),0.79,1.2,7e+3)
% quadgk(@(z)pfDESext2(z,x,u,J_jump,param),0.75,1.5)
