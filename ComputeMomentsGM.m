function [ muX,sigmaX,gammaX,kappaX ] = ComputeMomentsGM(mu,sigma,lambda)
%ComputeMomentsGM is a function for computing mean, variance, skweness and
%kurtosis of a univariate gaussian mixture random variable.
%   INPUT:
%      mu = vector of means
%      sigma = vector of standard deviations
%      lambda = vector of mixture proportions
%   OUTPUT:
%      muX = mean
%      sigma2X = variance
%      gammaX = skweness
%      kappaX = kurtosis

n = length(mu);
%% Mean
muX = 0;
for i = 1 : n
	muX = muX + lambda(i) * mu(i);
end

%% Variance
sigma2X = -muX^2;
for i = 1 : n
	sigma2X = sigma2X + lambda(i) * (mu(i)^2 + sigma(i)^2);
end
sigmaX = sqrt(sigma2X);
%% Skweness
gammaX = -3 * muX * sigma2X - muX^3;
for i = 1 : n
	gammaX = gammaX + lambda(i) * (mu(i)^3 + 3 * mu(i) * sigma(i)^2);
end
gammaX = gammaX / (sqrt(sigma2X))^3;

%% Kurtosis
kappaX = -muX^4 - 6 * muX^2 * sigma2X - 4 * gammaX * sqrt(sigma2X)^3 * muX;
for i = 1 : n
	kappaX = kappaX + lambda(i) * (mu(i)^4 + 6 * mu(i)^2 * sigma(i)^2 + 3 * ...
		sigma(i)^4);
end
kappaX = kappaX / sigma2X^2;
end

