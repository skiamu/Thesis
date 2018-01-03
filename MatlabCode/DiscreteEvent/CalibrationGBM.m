function [param] = CalibrationGBM(logReturns,dt)
%CalibrationGBM is a function for calibratin the two parameters of a
%Geometri Brownian Motion (GBM) on real data.
%   INPUT:
%      logReturns = market log-returns [vector]
%      dt = time step
%   OUTPUT:
%      param.mu =  drift estimation
%      param.sigma = volatility estimation

Sr = var(logReturns); % sample variance
r_bar = mean(logReturns); % sample mean
sigma_hat = sqrt(Sr / dt);
mu_hat = r_bar / dt + 0.5 * sigma_hat^2;
param.mu = mu_hat;
param.sigma = sigma_hat;
%% plot
% X = [0 cumsum((mu_hat - 0.5 * sigma_hat^2) * dt + sigma_hat * sqrt(dt) * randn([1 length(Returns)-1]))];
% figure
% plot(1:length(Returns),cumsum(Returns),'r',1:length(Returns),X,'b')

end % CalibrationGBM

