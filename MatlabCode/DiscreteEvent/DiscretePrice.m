function [D,DiscreteS] = DiscretePrice(S,J)
%DiscretePrice computes the discrite price of a given timeseries
%   INPUT:
%      S = price timeseries
%      J = jump size
%   OUTPUT:
%      D = vector of 1,0 or -1
%      DiscreteS = discrete price timeseries

N = length(S);
DiscreteS = S;
D = zeros([N 1]);
Baseline = S(1); % initial baseline is the initial price
for i = 2 : N
	r = S(i) / Baseline - 1; % return
	if abs(r) > J % abs return value greater than threshold
		Baseline = Baseline * (1 + sign(r) * J); % update baseline
		DiscreteS(i) = Baseline; % update discrete dynamics
		D(i) = sign(r);
	else
		D(i) = 0;
		DiscreteS(i) = DiscreteS(i-1); % discrete dynamics doesn't jump
	end
end
end % DiscretePrice


