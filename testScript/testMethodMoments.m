% test method of moment
% run script
clc; close all; clear variables;
addpath(genpath(pwd))
rng default
%% 1) get data from Yahoo
k = 2; M = 3;
Sample = [0.0324;0;0;3;0.0546;0.0445;-0.46;4.25;0.1062;0.1477;-0.34;5.51;0;0;0.0342];
SampleFreq = 'y';
[ x, error, lambda ] = MethodofMomentsGM( Sample,k,M,SampleFreq );
