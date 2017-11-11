function [param,CalibrationData] = GMcalibrationML( Returns,k,M )
%GMcalibrationML calibrate a gaussian mixture model using the maximum
%likelihood method
%   INPUT:
%      Returns = asset class returns [matrix]
%      k = number of mixture components [scalar]
%      M = asset allocation dimension [scalar]
%   OUTPUT:
%      param = struct array of parameters, each component of the array is a
%              struct with the following fields: mu, S, lambda
%      CalibrationData = struct with calibration informationeta = 0.01;
%   REMARKS :
%      1) problems for i = N (hence i = 1 : N -1)

eta = 0.01;
lambda = 0:eta:1;
N = length(lambda);
% the initial condition must be fealible (matrixes posotive definite)
d = rand(6);
x0 = [rand([6 1]); sqrt(diag(d'*d));-1 + (2).*rand([3 1])]; % smart initial point
% x0 = rand([15 1]); % naive initial point
x = zeros([15 N]);
f = zeros([N 1]);
lb = [-ones([k*M 1]); zeros([k*M 1]); -ones([3 1])]; % lower bound
ub = ones([15 1]); % upper bound
options = optimoptions(@fmincon,'Algorithm','interior-point',...
	'FiniteDifferenceType','central');
for i = 1 : N-1
	i
	[x(:,i),f(i)] = fmincon(@(x)-obj(x,Returns,k,lambda(i)),x0,[],[],[],[],lb,ub,...
		@(x)nonlinconst(x,k,lambda(i)),options);
% 	x0 = x(:,i);
end
f = -f;
[~,idxMax] = max(f);
x_star = x(:,idxMax);
param = MakeParam(x_star,k,lambda(idxMax));
CalibrationData.LogL = f(idxMax);
end


function f = obj(x,X,k,lambda)
% 1) compute LogL
param = MakeParam(x,k,lambda); % convert vector x into the struct of parameters
try
f = sum(log(GMdensity(X,param,k))); % compute LogL
catch % matrix not positive definite
	f = -1e7;
end
end

function [c, ceq] = nonlinconst(x,k,lambda)
%nonlinconst set the positive-definitness and unimodality constraint
%   INPUT:
%      x = 
%      k = 
%      lambda = 
%   OUTPUT:
%      c = 
%      ceq = 

param = MakeParam(x,k,lambda);
% 1) positive-definite
p = zeros([k 1]);
for i = 1 : k
	[~, p(i)] = chol(param(i).S);
end
% 2) unimodality
q = zeros([3 1]);
for i = 1 : 3
	q(i) = (param(2).mu(i) - param(1).mu(i))^2  - 27 * param(1).S(i,i) * ...
		param(2).S(i,i) / (4*(param(1).S(i,i) + param(2).S(i,i)));
end
% 3) set the constraints
ceq = p;
c = q;
end

function param = MakeParam(x,k,lambda)
mu = [x(1:3) x(4:6)];
sigma = [x(7:9) x(10:12)];
Rho = x(13:15);
lambda = [lambda; 1 - lambda];
RhoM = [1 Rho(1) Rho(2); Rho(1) 1 Rho(3); Rho(2) Rho(3) 1];
param(k) = struct();
for i = 1 : k
	param(i).mu = mu(:,i);
	param(i).S = corr2cov(sigma(:,i),RhoM);
	param(i).lambda = lambda(i);
end
end