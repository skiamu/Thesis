function [Binomial, tau] = SimulationED(param,Nsim, Nstep,J_jump,model)
%SimulationED simulated the random variables ih the event-driven dynamics
%for the basic model and extension1
%   INPUT:
%      param = model parameters [struct]
%      Nsim = number of MC simulation [scalar]
%      Nstep = number of time steps [scalar]
%      J_jump = jump size [scalar]
%      model = {'basic','ext1'}
%  OUTPUT:
%     Binomial = jump sign [matrix Nsim x Nstep]
%     tau = holding times [matrix Nsim x Nstep]

p = param.p;
%% simulation tau
if strcmp(model,'basic')
	lambda = param.lambda;
	tau = exprnd(1 / lambda, [Nsim Nstep]);
else % model ext1
	mu = param.mu; sigma = param.sigma;
	mu_tilde = mu-sigma^2/2;
	F_tau = @(t)1-(2 * cosh(mu_tilde * J_jump / sigma^2) * Kappa(t,mu_tilde,sigma,J_jump));
	t = (1e-4:1e-5:4)';
	[u,idxUnique] = unique(F_tau(t)); % some value may repete themselves. interp1 raises an error
	u_tilde = rand([Nsim Nstep]);
	tau = interp1(u,t(idxUnique),u_tilde);
end

%% simulation Binomial
Binomial = -ones([Nsim Nstep]);
Binomial(rand([Nsim Nstep]) < p) = 1;
%Binomial = Binomial(randperm(Nsim),randperm(Nstep));
end % SimulationED

function [f] = Kappa(y,mu_tilde,sigma,J_jump)
%   INPUT:
%      y = (z+-csi) / x (column vector)
%      x = ptf value
%      mu_tilde = mu - 0.5 * sigma^2;
%      sigma = volatility
%      r = cash yearly return
%      J_jump = jump size

epsilon = 1e-8; % series truncation tolerance
n = length(y);
t = min(y);
NN = sqrt(max(1, -8 * J_jump^2 ./ (pi^2 * sigma^2 * t) .* ...
	(log((pi^3 * sigma^2 * t * epsilon) / (16 * J_jump^2))...
	- J_jump * mu_tilde / sigma^2))); % number of series terms

global freeMemory
if 8 * NN * n * 1e-6 >= freeMemory % check if there's enogh memory
	f = [0; Kappa(y(2:end),mu_tilde,sigma,J_jump)];
else
	N = 1:NN;
	% y is a column vector, N must be a row vector. In this case the operation
	% y.^N gives a matrix [lenght(y) length(N)]. To get the sum of the partial
	% sum we need to sum by rows (e.g. sum(,2))
	f = sigma^2 * pi / (4 * J_jump^2) * (sum(N .* (-1).^(N+1) ./ ...
		(mu_tilde^2 / (2 * sigma^2) + (sigma^2 * N.^2 * pi^2) / (8 * J_jump^2)) .*  ...
		exp(-(mu_tilde^2 / (2 * sigma^2) + (sigma^2 * N.^2 * pi^2) / (8 * J_jump^2)) .* y)...
		.* sin(pi * N / 2),2));
end

end % Gamma



