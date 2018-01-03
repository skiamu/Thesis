% test script for Extensions
clc; close all; clear variables;
addpath(genpath(pwd))
load('stocks.mat')
%% discrete dymanics
freq = 'd';
M = 3;
[Returns,~,stocks] = getReturns( freq, M );
LogReturns = log(1 + Returns(:,3));
dt = 1 / 252; % time expressed in years
[mu_hat,sigma_hat] = CalibrationGBM(LogReturns,dt,1);