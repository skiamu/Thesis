clc; close all; clear variables;
RhoM = [1 0.0633 0.0207;
	0.0633 1 -0.0236;
	0.0207 -0.0236 1];
param(1).mu = [0.000611 0.001373 0.002340]';
param(2).mu = [0.000683 -0.016109 -0.017507]';
param(1).S = corr2cov([0.000069 0.00566 0.019121],RhoM);
param(2).S = corr2cov([0.000062 0.006168 0.052513],RhoM);
param(1).lambda = 0.98;
param(2).lambda = 0.02;
n = length(param);
W = 0;
for i = 1 : n
	W = W + param(i).lambda * param(i).S;
	for j = 1 : i-1
		W = W + param(i).lambda * param(j).lambda * ...
			(param(i).mu - param(j).mu) * (param(i).mu - param(j).mu)';
	end
end

u1 = [0.215403049464726;0.0521792032793122;0.732417747255962];
u2 = round([0;0.274255267162805;0.725744732837195],3);
u3 = [0;0.274255267162805;0.725744732837195];
u4 = [0.266900956843266;0;0.733099043156735];
s2 = sqrt(u2' * W * u2);
s3 = sqrt(u3' * W * u3);
s4 = sqrt(u4' * W * u4);
s4 > s3
s4 > s2
alpha = 0.01;
VaR = 0.07;
sigmaMax2 = (VaR / norminv(1-alpha))^2;
% VaRGM(param,u,alpha);
% quadgk(@(alpha)VaRGM(param,u,alpha),0,alpha) / (alpha)

mu = [param(1).mu(2);param(2).mu(2)];
sigma = sqrt([param(1).S(2,2);param(2).S(2,2)]);
lambda = [param(1).lambda,param(2).lambda]
[muX,sigmaX,gammaX,kappaX ] = ComputeMomentsGM(mu,sigma,lambda)
sigmaX^2








