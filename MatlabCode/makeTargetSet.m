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
X{end} = ((1+theta)^1:eta:UB)';
end



