function [ x, error, lambda, param ] = MethodofMomentsGM( Sample,k,M,SampleFreq )
%MethodofMomentsGM is a function for calibrating a Gaussian Mixture models
%by moment metching. It is used when time-series are not available.
%   INPUT:
%      Sample = vector of sample moments. Sample = [muC,sigmaC,gammaC,kappaC,...
%               ...,muE,sigmaE,gammaE,kappaE, rho12,rho13,rho23].
%
%      k = number of gaussian components
%      M = dimension
%   OUTPUT:
%      x =
%      error =
%   REMARKS: data in sample have an annual frequancy, they must be
%   converted to weekly

if strcmp(SampleFreq,'y')
	Sample([1 5 9]) = (1 + Sample([1 5 9])).^(1/52) - 1;
	Sample([2 6 10]) = Sample([2 6 10]) / sqrt(52);
end
eta = 0.01;
lambda = 0:eta:1;
error = zeros(length(lambda),1);
x = zeros(length(lambda),15);
lb = [-ones([k*M 1]); zeros([k*M 1]); -ones([3 1])];
ub = [ones([k*M 1]);ones([k*M 1]);ones([3 1])];
% x0 = [0.000611; 0.001373;0.002340;0.000683;-0.016109;-0.017507;0.000069;...
% 	0.00566;0.019121;0.000062;0.006168;0.052513;0.0633;0.0207;-0.0236];
x0 = rand([15 1]);
for i = 1 : length(lambda)
	i
	[x(i,:), error(i)] = lsqnonlin(@(x) obj(x,Sample,k,M,lambda(i)),x0,lb,ub);
% 	x0 = x(i,:);
end


if nargout > 3
	[~,idxMin] = min(error);
	lambda = lambda(idxMin);
	x = x(idxMin,:);
	disp(['calibration error = ', num2str(error(idxMin))]);
	param = cell([k 1]);
	mu1 = x(1:3); mu2 = x(4:6);
	sigma1 = x(7:9); sigma2 = x(10:12);
	Rho = x(13:15);
	RhoM = [1 Rho(1) Rho(2); Rho(1) 1 Rho(3); Rho(2) Rho(3) 1];
	S1 = corr2cov(sigma1,RhoM);
	S2 = corr2cov(sigma2,RhoM);
	param{1} = {mu1',S1,lambda}; % traspose since column vector are wanted
	param{2} = {mu2',S2,1-lambda};
end

end %MethodofMomentsGM

function F = obj(x,Sample,k,M,lambda)

% 1) extract parameters
% Mean = [mu1; mu2]
Mean = x(1:k*M); Mean = reshape(Mean, [M k])';
% Std = [std1; std2]
Std = x(k*M+1:2*k*M) ; Std = reshape(Std, [M k])';
% Rho = [rho12, rho13, rho23]
Rho = x(2*k*M+1:end);
lambda = [lambda; 1-lambda];
F = zeros(size(x)); % initialization

% 2) computes theoretical moments
for i = 1 : M
	mu = Mean(:,i); % first component's means
	sigma = Std(:,i); % first component's variances
	[muX,sigmaX,gammaX,kappaX] = ComputeMomentsGM(mu,sigma,lambda);
	F(1+(i-1)*4:i*4) = [muX;sigmaX;gammaX;kappaX] - Sample(1+(i-1)*4:i*4);
end

% 3) compute correlations
RhoSample = Sample(13:end);
RhoM = [1 RhoSample(1) RhoSample(2); RhoSample(1) 1 RhoSample(3); RhoSample(2) RhoSample(3) 1];
Sigma = corr2cov(Sample([2 5 9]),RhoM);

k = 13;
for i = 1 : 3
	for j = 1 : i-1
		F(k) = lambda(1) * Rho(k-12) * Std(1,j) * Std(1,i) + lambda(2) * Rho(k-12) * Std(2,j) * Std(2,i) + ...
			lambda(1) * lambda(2) * (Mean(1,i) - Mean(2,i)) * (Mean(1,j) - Mean(2,j)) - Sigma(i,j);
		k = k + 1;
	end
end
% F = SetConstraints(F,Mean,Std,Rho);

end

function F = SetConstraints(F,Mean,Std,Rho)
% 1) check for unimodality for each marginal
Unimodality = zeros([3 1]);
for i = 1 : 3
	Unimodality(i) = (Mean(1,i)-Mean(2,i))^2 < 27 * Std(1,i)^2 * ...
		Std(2,i)^2 / (4 * (Std(1,i)^2 + Std(2,i)^2));
end

% 2) check for positive-definitness
RhoM = [1 Rho(1) Rho(2); Rho(1) 1 Rho(3); Rho(2) Rho(3) 1];
S1 = corr2cov(Std(1,:),RhoM);
S2 = corr2cov(Std(2,:),RhoM);
[~,p1] = chol(S1);
[~,p2] = chol(S2);

if any(~Unimodality) || any([p1 p2])
	F = 1e10 .* F;
end

end










