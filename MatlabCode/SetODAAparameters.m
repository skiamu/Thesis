% this file set the parameters for the ODAA algorithm
VaR = 0.07; % monthly
alpha = 0.01; % confidence level VaR
switch freq % number of time step for a 2-year investment
	case 'wk'
		N = 104;
		NstepPlot = 26;
		VaR = VaR / 2;
	case 'm'
		N = 24;
		NstepPlot = 6;
	case 'q'
		N = 8;
		NstepPlot = 3;
		VaR = VaR * sqrt(3);
end
theta = 0.07; % yearly target return
eta = 1e-3; % target set discretization
[ X ] = makeTargetSet(N,theta,eta);


function [ X ] = makeTargetSet(N,theta,eta)
%makeTargetSet creates the discretized target sets used in the DPalgorithm
%   INPUT:
%      N = nuber of time step
%      theta = yearly target return
%      eta = discretization step
%   OUTPUT:
%      X = cell array of target sets

X = cell([N+1 1]); % initialization
LB = 0.5; % lober bound approximation
UB = 1.9; % upper bound approximation
for i = 2 : N
	X{i} = (LB:eta:UB)';
end
X{1} = 1; 
X{end} = ((1+theta)^2:eta:UB)';
end