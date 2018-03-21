% run script Event-driven approach

%% 1) read data
clc; close all; clear variables;
global freeMemory 
freeMemory = 0.95 * get_free_mem();
addpath(genpath(pwd))
DefaultPlotting
load('LIBOR.mat') % risk-free asset
LIBOR = LIBOR(1:2000) / 100;
load('SP500Fut.mat') % risky asset
idx1 = find(DateSP500Fut == '22-Jan-2010');
idx2 = find(DateSP500Fut == '25-Apr-2016');
S = SP500Fut(idx1:idx2);
% S = SP500Fut;
Returns = S(2:end) ./ S(1:end-1) - 1; % compute risky asset returns
save('/home/andrea/Thesis/PythonCode/Event_Driven/Data.mat','S','Returns','LIBOR')
%%  2) Calibration
model = 'ext1'; % select from {'basic','ext1','ext2'}
J_jump = 0.1; % jump size
dt = 1 / 252; % time in years
switch model
	case 'basic'
		[param] = Calibration(S,model,J_jump,dt);
		p = param.p; lambda = param.lambda; param.r = 0.03;
	case 'ext1'
		[param] = Calibration(Returns,model,J_jump,dt);
		param.r = 0.01; % cash return
		mu = param.mu; sigma = param.sigma;
		mu_tilde = mu - 0.5 * sigma^2;
		param.p = (exp(2 * mu_tilde * J_jump / sigma^2) - 1 ) / ...
			(2 * sinh(2 * mu_tilde * J_jump / sigma^2)); % probability positive jump
	case 'ext2'
		MarketData.LIBOR;
		MarketData.RiskyAsset = S;
		[param] = Calibration(MarketData,model,J_jump,dt);
end

%% 3) ODAA algorithm
SetODAAParamED % set parameter for the ODAA algorithm
tic
[U,J] = ODAAalgorithmDES(N,X,J_jump,param,model);
time = toc
p_star = J{1};

%% 4) Validation
Nsim = 1e+6;
[p_starMC,Times] = ValidationDES(X,U,Nsim,N,param,J_jump,model);

%% 5) print the maps
% plot the maps
i=6;
for k = 10 : -2 : 2
	idx = find(X{k} <= 1.9 & X{k} >= 0.6);
	subplot(3,2,i)
	area(X{k}(idx),U{k}(idx))
	xlim([0.8 1.9]);
	xlabel('portfolio value')
	title(strcat('k = ',num2str(k-1)))
	grid on
	i = i - 1;
end

% print the maps and save the results
print(strcat('/home/andrea/Thesis/Latex/final/Images/maps',model),'-dpng', '-r900');
save(strcat('/home/andrea/Thesis/MatlabCode/DiscreteEvent/',model,'.mat'),'J',...
	'U', 'X', 'p_star','p_starMC','time','Times')

%% 6) check the density
u = 0.5;
x = 0.7;
csi = x * J_jump * u;
z = (0.5:1e-5:1.2)';
f = pfDES(z,x,u,J_jump,param);
plot(z,f,'r.-')
xlim([0.6 1.9]);
tic
1-trapz(z,pfDES(z,x,u,J_jump,param))
toc

% figure
% for u = -1:0.5:1
% 	plot(z,pfDESext1(z,x,u,J_jump,param),'.-','DisplayName',['u = ', num2str(u)])
% 	hold on
% end
% legend('show')
% title('Portfolio value density')
% print('/home/andrea/Thesis/Latex/final/Images/PtfDensity','-dpng', '-r900');

