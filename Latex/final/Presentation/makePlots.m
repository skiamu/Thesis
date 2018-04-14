%% Raggiungibilità stocastica
mu = 0.1;
sigma = 0.30;
dt = 1/252;
Nsim = 252*4;
DeltaW = randn([1 Nsim]);
S = exp(mu*dt+sigma*sqrt(dt)*cumsum(DeltaW));
plot(S)


%% scatterplot mistura gaussiana
mu1 = [2 5];          % Mean of the 1st component
sigma1 = [1 0; 0 3]; % Covariance of the 1st component
mu2 = [-1 -2];        % Mean of the 2nd component
sigma2 = [1 0; 0 3];  % Covariance of the 2nd component
rng('default') % For reproducibility
r1 = mvnrnd(mu1,sigma1,10000);
r2 = mvnrnd(mu2,sigma2,1700);
X = [r1; r2];
gm = fitgmdist(X,2)
scatter(X(:,1),X(:,2),10,'.') % Scatter plot with points of size 10
hold on
ezcontour(@(x,y)pdf(gm,[x y]),[-5 6 -8 13])
title('')
print('/home/andrea/Thesis/Latex/final/Presentation/Images/Scatter','-dpng', '-r900');

%% Exit time

mu = 0.07;
sigma = 0.16;
dt = 1/252;
Nsim = ceil(252*0.3);
DeltaW = randn([1 Nsim]);
X = [0 mu*dt+sigma*sqrt(dt)*cumsum(DeltaW)];
plot((0:Nsim)/252,X)
grid on


print('/home/andrea/Thesis/Latex/final/Presentation/Images/ExitTime','-dpng', '-r900');


