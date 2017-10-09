% test MCECMalgorithm
clc;clear variables; close all;
rng(1999)
load('stocks.mat') % load stocks data
load('dist.mat') % load parameters computed in R
load('Yahoo.mat')
% X = smi_stocks(:,1:2);
X = Returns;
[N, d] = size(X);
toll = 1e-4; 
maxiter = 500;
GHmodel = 'NIG';
tic
[theta,LogL,exitFlag,numIter] = MCECMalgorithm(toll,maxiter,X,GHmodel);
time = toc;

theta{numIter}

lambda = theta{numIter}{1}; Chi = theta{numIter}{2};
Psi = theta{numIter}{3}; mu = theta{numIter}{4}; Sigma = theta{numIter}{5};
gamma = theta{numIter}{6};
alpha_bar = sqrt(Psi * Chi)
