% set ODAA parameters in the event-driven case

N = 10; % number of events
theta = 0.07; % yearly target return
eta = 1e-3/5; % target set discretization
n = 3;
[ X ] = makeTargetSet(N,theta,eta,n);
function [ X ] = makeTargetSet(N,theta,eta,n)
%makeTargetSet creates the discretized target sets used in the DPalgorithm
%   INPUT:
%      N = nuber of time step
%      theta = yearly target return
%      eta = discretization step
%   OUTPUT:
%      X = cell array of target sets

X = cell([N+1 1]); % initialization
LB = 0.8; % lober bound approximation
UB = 2; % upper bound approximation
for i = 2 : N
	X{i} = (LB:eta:UB)';
end
X{1} = 1; 
X{end} = ((1+theta)^n:eta:UB)';
end