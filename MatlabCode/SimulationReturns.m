function [ w ] = SimulationReturns(param,Nsim,M,Nstep,model)
%SimulationReturns is a function for simulating asset class returns
%according to the model specified in input
%   INPUT:
%      param = cell array or struct of model parameters
%      Nsim = number of MC simulations
%      M = asset class dimension
%      Nstep = number of time intervals
%      SimulationMethod =
%      model = 'Gaussian', 'Mixture', 'GH'
%   OUTPUT:
%      w =
% REMARK : gestire l'input GM
switch model
	case 'Gaussian'
		[ w ] = SimulationG(param,Nsim,Nstep,M);
	case 'Mixture'
		[ w ] = SimulationGM(param,Nsim,Nstep,M);
	case {'GH','NIG','t'}
		[ w ] = SimulationGH(param,Nsim,Nstep,M);
	otherwise
		error('invalid model %s',model);
end

end

function [ w ] = SimulationG(param,Nsim,Nstep,M)
mu = param.mu; Sigma = param.S;
w = zeros([Nsim M Nstep]); %initialization
for k = 1 : Nstep
	w(:,:,k) = mvnrnd(mu,Sigma,Nsim);
end
end % SimulationG

function [ w ] = SimulationGM(param,Nsim,Nstep,M)
%SimulationGM is a function for simulating asset class returns over a
%finite time grid
%   INPUT:
%      GM = fitgmdist's output, it contains all the information about the
%      fitted gaussian mixture distribution
%      param = cell array with distrubution's paraeters
%      Nsim = numer of MC simulation
%      M = asset allocation dimension
%      Nstep = number of time intervals
%      simulationMethod = 'built-in', 'Normal'
%   OUTPUT:
%      w =

w = zeros([Nsim M Nstep]); %initialization
for k = 1 : Nstep
	X = [];
	for i = 1 : length(param)
		X = [X; mvnrnd(param(i).mu,param(i).S, round(Nsim * param(i).lambda))];
	end
	X = X(randperm(end),:);
	w(:,:,k) = X;
end
end
% SimulationGM

function [ w ] = SimulationGH(param,Nsim,Nstep,M)
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
