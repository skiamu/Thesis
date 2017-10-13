function [ w ] = SimulationGH(param,Nsim,M,Nstep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% extract parameters
lambda = param.lambda; Chi = param.Chi; Psi = param.Psi;
mu = param.mu; Sigma = param.sigma; gamma = param.gamma;

w = zeros([Nsim M Nstep]); %initialization
Z = randn([Nsim M Nstep]); % normal gaussian

for k = 1 : Nstep
	W = gigrnd(lambda, Psi, Chi, Nsim); % simulate GIG
	W = W *  ones([1 M]); % triplicate the column
	w(:,:,k) = repmat(mu',[Nsim 1]) + W .* repmat(gamma',[Nsim 1]) + sqrt(W) .* ...
		(Z(:,:,k) * chol(Sigma));
end

end % SimulationGH

