function [param] = CalibrationBasic(AssetPrice,J_jump,dt)
%CalibrationBasic id a function for calibrating the basic model (parameters p and lambda)
%   INPUT:
%      AssetPrice = risky-asset prices [vector]
%      J_jump = jump size
%   OUTPUT:
%      param = struct of model parameter

%% 1) discrete dymanics
[D,DiscreteS] = DiscretePrice(AssetPrice,J_jump);
figure 
xx = (1:length(D))';
plot(xx,AssetPrice,'b-',xx,DiscreteS,'r-')

%% 2) calibration
LB = [0;0];
UB = [1;Inf];
X0 = [rand([100 1]) rand([100 1])* 20 ]; % initialize different initial points
x_star = zeros([100 2]); f_star = zeros([100 1]);
for i = 1 : 100
	[x_star(i,:), csi] = fmincon(@(x)-obj(D,x,dt),X0(i,:),[],[],[],[],LB,UB,[],[]);
	f_star(i) = -csi;
end
[~,idxMax] = max(f_star);
x_star = x_star(idxMax,:);
param.p = x_star(1);
param.lambda = x_star(2);

end % CalibrationBasic

function f = obj(D,x,dt)
%obj define the log-likelihood function to be maximized
%   INPUT:
%      D = jump data
%      x = vector of parameters [p,lambda]
%      dt = time step
%   OUTPUT:
%      f = log-likelihood function

p = x(1); lambda = x(2); % extract parameters
Nzero = sum(D == 0);
Nuno = sum(D == 1);
NmenoUno = sum(D == -1);
assert(Nzero + Nuno + NmenoUno == length(D));
f = -Nzero * lambda*dt + Nuno * log((1-exp(-lambda * dt)) * p) +...
	NmenoUno * log((1-exp(-lambda*dt)) * (1-p)); % log-likelihood function
end % obj
