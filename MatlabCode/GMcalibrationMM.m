function [param,CalibrationData] = GMcalibrationMM(Returns,k,M )
%GMcalibrationMM calibrates a gaussian mixture model by using the method of
%moments
%   INPUT:
%      Returns = asset class returns [matrix]
%      k = number of mixture components
%      M = asset allocation dimension
%   OUTPUT:
%      param = struct array of parameters, each component of the array is a
%              struct with the following fields: mu, S, lambda
%      CalibrationData = struct with calibration information
%   REMARKS:
%      1) this function is not general. It works only with mixtures of 2
%      components and in dimension 3
%      2) by using a rolling initial condition in lsqnonlin we get better results. Try
%      to use both methods and see the differences in term of allocation
%      maps
%% 1) compute sample statistics
Sample = zeros([4 M]); % each column will contain mean, std, skwe and kurtosis
Sample(1,:) = mean(Returns);
Sample(2,:) = std(Returns);
Sample(3,:) = skewness(Returns);
Sample(4,:) = kurtosis(Returns);
Corr = corr(Returns); % sample correlation matrix
SampleCorr = [Corr(1,2); Corr(1,3); Corr(2,3)];
clear Corr;

debug = false; % set to true 
if debug
	Sample = [0.0324;0.00;0;3;
		0.0546;0.0445;-0.46;4.25;
		0.1062;0.1477;-0.34;5.51];
	Sample([1 5 9]) = (1 + Sample([1 5 9])).^(1/52) - 1;
	Sample([2 6 10]) = Sample([2 6 10]) / sqrt(52);
	Sample = reshape(Sample, [4 3]);
	SampleCorr = [0;0;0.0342];
end
%% 2) least squares optimization
eta = 0.01;
lambda = 0:eta:1;
N = length(lambda);
error = zeros([N 1]); % initialization
x = zeros([15 N]); % initialization
lb = [-ones([k*M 1]); zeros([k*M 1]); -ones([3 1])]; % lower bound
ub = ones([15 1]); % upper bound
R = corr(Returns);
x0 = [repmat(mean(Returns)',[2 1]) + (1e-2 * randn([6 1]));
	repmat(sqrt(diag(cov(Returns))),[2 1]) + (1e-2 * randn([6 1]));
	R(1,2);R(1,3);R(2,3)];
clear R;
% d = rand(6);
% x0 = [rand([6 1]); sqrt(diag(d'*d));-1 + (2).*rand([3 1])]; % smart initial point
% x0 = rand([15 1]);
options = optimoptions(@lsqnonlin);
for i = 1 : N
	i
	[x(:,i), error(i)] = lsqnonlin(@(x) obj(x,Sample,SampleCorr,k,M,lambda(i)),...
		x0,lb,ub,options);
% 	x0 = x(:,i); % rolling initial point
end
%% 3) choose best lambda
% we select tha lambda which minimizes the residual error
[minError,idxMin] = min(error);
x_star = x(:,idxMin); % vector of optimal parameters
lambda = lambda(idxMin); 
CalibrationData.error = minError; % return calibration information

%% 4) assemble the struct param
param = MakeParam(x_star,k,lambda);

%% 5) compute LogL
CalibrationData.LogL = sum(log(GMdensity(Returns,param,k)));

end % GMcalibrationMM


function F = obj(x,Sample,SampleCorr,k,M,lambda)
%obj if the objective function for the least squares problem, it's a system
%of moment equations.
%   INPUT:
%      x = vector of parameters(1:6 means 7:12 stdev 13:15 correlations)
%      Sample = matrix of sample moments (each colums contains mean,std,skw and kurtosis)
%      SampleCorr = vector of sample cross correlation (rho12,rho13,rho23)
%      k = number of gaussian components
%      M = asset allocation dimension
%      lambda = proportion
%   OUTPUT:
%      F = vector objective function

%% 1) extract parameters
% for readability purposes with extract the parameters from the column
% vector x. Mean is a matrix whose columns are mean vectors of gaussian
% component and Std collects standard deviation of each gaussian component
Mean = reshape(x(1:k*M), [M k]);
Std = reshape(x(k*M+1:2*k*M), [M k]);
Rho = x(2*k*M+1:end);
lambda = [lambda; 1-lambda];

%% 2) compute theoretical moments
F = [];
for i = 1 : M
	marginalMean = Mean(i,:); % means of the i-th marginal
	marginalStd = Std(i,:); % std of the i-th marginal
	[muX,sigmaX,gammaX,kappaX] = ComputeMomentsGM(marginalMean,marginalStd,lambda);
	F = [F; [muX;sigmaX;gammaX;kappaX] - Sample(:,i)];
end

%% 3) compute theoretical correlations
SampleCorrMatrix = [1 SampleCorr(1) SampleCorr(2);
	SampleCorr(1) 1 SampleCorr(3);
	SampleCorr(2) SampleCorr(3) 1];
SampleCovMatrix = corr2cov(Sample(2,:),SampleCorrMatrix);
k = 13;
for i = 1 : 3
	for j = 1 : i-1
		F(k) = lambda(1)*Rho(k-12)*Std(i,1)*Std(j,1) + ...
			lambda(2)*Rho(k-12)*Std(i,2)*Std(j,2) + ...
			lambda(1)*lambda(2)*(Mean(i,1)-Mean(i,2))*(Mean(j,1)-Mean(j,2)) -...
			SampleCovMatrix(i,j);
		k = k + 1;
	end
end
%% 4) check for positive-defitness and unimodality
F = SetConstraints(F,Std,Rho,x,lambda);

end % obj


function F = SetConstraints(F,Std,Rho,x,lambda)
RhoM = [1 Rho(1) Rho(2); Rho(1) 1 Rho(3); Rho(2) Rho(3) 1];
S1 = corr2cov(Std(:,1),RhoM);
S2 = corr2cov(Std(:,2),RhoM);
[~,p1] = chol(S1);
[~,p2] = chol(S2);

param = MakeParam(x,2,lambda);
unimodality = zeros([3 1]);
for i = 1 : 3
	unimodality(i) = (param(2).mu(i) - param(1).mu(i))^2  < 27 * param(1).S(i,i) * ...
		param(2).S(i,i) / (4*(param(1).S(i,i) + param(2).S(i,i)));
%    q(i) = abs(param(1).mu(i)-param(2).mu(i)) - min(param(1).S(i,i),param(2).S(i,i));
end

if any([p1 p2] ~= 0) || any(~unimodality)
	F = F .* 1e+10;
end
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











