function [ X ] = makeTargetSet(N,theta,eta )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

X = cell([N+1 1]); % initialization
LB = 0.5;
UB = 1.9;
for i = 2 : N
	X{i} = (LB:eta:UB)';
end
X{1} = 1; X{end} = ((1+theta)^1:eta:UB)';

end

