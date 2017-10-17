% test MCECMalgorithm
clc;clear variables; close all;
rng(1999)
addpath(genpath(pwd))
load('stocks.mat') % load stocks data
load('Yahoo.mat')
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

%% new functions
% lambda = -0.5; Chi = 0.7; Psi = 0.7;
% [mu, lambda, gamma, Sigma, chi, psi]=multiGH_mcecm_clam_fit(X,lambda,Chi,Psi,maxiter,toll);
