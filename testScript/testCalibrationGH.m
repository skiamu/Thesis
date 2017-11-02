% test MCECMalgorithm
clc;clear variables; close all;
rng(1999)
addpath(genpath(pwd))
load('stocks.mat') % load stocks data
load('Yahoo.mat')
X = Returns;
[N, d] = size(X);
toll = 1e-10; 
maxiter = 2000;
GHmodel = 'VG';
tic
[theta,LogL,exitFlag,numIter] = MCECMalgorithm_VG(toll,maxiter,X,GHmodel);
time = toc;

theta{numIter}
% sqrt(theta{numIter}{2}*theta{numIter}{3})
lambda = theta{numIter}{1};
alpha = theta{numIter}{2};
mu = theta{numIter}{3};
Sigma = theta{numIter}{4};
gamma = theta{numIter}{5};
Chi = theta{numIter}{6};
Psi = theta{numIter}{7};

%% new functions
% lambda = -0.5; Chi = 0.7; Psi = 0.7;
% [mu, lambda, gamma, Sigma, chi, psi]=multiGH_mcecm_clam_fit(X,lambda,Chi,Psi,maxiter,toll);
