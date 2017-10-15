function [x] = MethodMomentsGM(X)
%MethodMomentsGM is a function for the calibration of a Gaussian Mixture
%(GM) model by using the method of moments.
%   INPUT:
%      X =  data matrix (asset class returns)
%   OUTPUT:
%      x = vector of parameters [lambda, mu1_1, mu2_1, mu3_1, sigma1_1,
%          sigma2_1,sigma3_1, mu1_2, mu2_2, mu3_2,sigma1_2, sigma2_2, sigma3_2,
%          rho12, rho13, rho23]

[~, d] = size(X);
% 1) compute sample moments
C = cov(X); % correlation matrix
MM = zeros([5 d]);
for i = 1 : 5
	MM(i,:) = mean(X.^i);
end

% 2) solve non-linear system
lb = [0, -1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, 0, -1, -1, -1];
ub = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Niter = 50;
x0 = 0.1 * randn([16 Niter]);
% x = zeros([16 Niter]);
% options = optimoptions(@fmincon, 'Algorithm','interior-point');
% for i = 1 : Niter
% 	i
% % 	x(:,i) = fmincon(@(x)0,x0(:,i),[],[],[],[],lb,ub,@(x)fminconstr(x,MM,C),options);
% % 	x(:,i) = fsolve(@(x)fbnd(x,MM,C),x0(:,1));
%    [x(:,i),r] = lsqnonlin(@(x)fbnd(x,MM,C),x0(:,i),lb,ub);
% 	r
% end
[x] = lsqnonlin(@(x)fbnd(x,MM,C),x0(:,1),lb,ub);
end % MethodMomentsGM

function [c,ceq] = fminconstr(x,MM,C)
c = []; % no nonlinear inequality
ceq = fbnd(x,MM,C); % the fsolve objective is fmincon constraints
end % fminconstr


function F = fbnd(x,MM,C)
%fbnd sets moment eqiations as constraints of the optimization problem

% 1) extract parameters for readability reasons
lambda = x(1);
mu1_1 = x(2);  mu2_1 = x(3);  mu3_1 = x(4);
sigma1_1 = x(5);  sigma2_1 = x(6);  sigma3_1 = x(7);
mu1_2 = x(8);  mu2_2 = x(9);  mu3_2 = x(10);
sigma1_2 = x(11);  sigma2_2 = x(12);  sigma3_2 = x(13);
rho12 = x(14);  rho13 = x(15); rho23 = x(16);

% 2) set moment equation
k = 0;
F = zeros([16 1]);
for i = 1 : 3 % for each asset class
	k = k + 1;
	mu1 = x(1+i); sigma1 = x(4+i);
	mu2 = x(7+i); sigma2 = x(10+i);
	F(k) = lambda * (mu1) + (1-lambda) * (mu2) - MM(1,i);
	F(k+1) = lambda * (sigma1^2 + mu1^2) + (1-lambda) * (sigma2^2 + mu2^2)  - MM(2,i);
	F(k+2) = lambda * (mu1^3 + 3 * mu1 * sigma1^2) + ...
		(1-lambda) * (mu2^3 + 3 * mu2 * sigma2^2) - MM(3,i);
	F(k+3) = lambda * (mu1^4 + 6 * mu1^2 * sigma1^2 + 3 * sigma1^4) + ...
		(1-lambda) * (mu2^4 + 6 * mu2^2 * sigma2^2 + 3 * sigma2^4) - MM(4,i);
end

% 3) set correlation equations
F(13) = lambda * rho12 * sigma1_1 * sigma2_1 + (1-lambda) * rho12 * sigma1_2 * sigma2_2 + ...
	lambda * (1-lambda) * (mu1_1 - mu1_2) * (mu2_1 - mu2_2) - C(1,2);

F(14) = lambda * rho13 * sigma1_1 * sigma3_1 + (1-lambda) * rho13 * sigma1_2 * sigma3_2 + ...
	lambda * (1-lambda) * (mu1_1 - mu1_2) * (mu3_1 - mu3_2) - C(1,3);

F(15) = lambda * rho23 * sigma2_1 * sigma3_1 + (1-lambda) * rho23 * sigma2_2 * sigma3_2 + ...
	lambda * (1-lambda) * (mu2_1 - mu2_2) * (mu3_1 - mu3_2) - C(2,3);

% 4) set last equation

F(16) = lambda * (mu1_1^5 + 10 * mu1_1^3 * sigma1_1^2 + 15 * mu1_1 * sigma1_1^4) + ...
	(1 - lambda) * (mu1_2^5 + 10 * mu1_2^3 * sigma1_2^2 + 15 * mu1_2 * sigma1_2^4) - MM(5,1);
end % fbnd









