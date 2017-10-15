function [ w ] = SimulationReturns(param,Nsim,M,Nstep,model,simulationMethod,GM )
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
		[ w ] = SimulationGM(GM,param,Nsim,M,Nstep,simulationMethod);
	case 'GH'
		[ w ] = SimulationGH(param,Nsim,M,Nstep);
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

function [ w ] = SimulationGM(GM,param,Nsim,M,Nstep,simulationMethod  )
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
switch simulationMethod
	case 'built-in'
		w = zeros([Nsim M Nstep]); %initialization
		for k = 1 : Nstep
			w(:,:,k) = random(GM,Nsim);
		end
	case 'Normal'
		% REMARK: this implementation work only if there are 2 mixture
		% components.
		
		% extract distribution parameters
		p1 = param{1}{3}; mu1 = param{1}{1}; Sigma1 = param{1}{2};
		mu2 = param{2}{1}; Sigma2 = param{2}{2};
		X1 = zeros([Nsim M Nstep]); X2 = zeros([Nsim M Nstep]);
		for k = 1 : Nstep
			X1(:,:,k) = mvnrnd(mu1,Sigma1,Nsim);
			X2(:,:,k) = mvnrnd(mu2,Sigma2,Nsim);
		end
		U = rand([Nsim M Nstep]);
		w = X2;
		idxFromX1 = find(U <= p1); % LINEAR indexes elements from X1
		w(idxFromX1) = X1(idxFromX1);
	otherwise
		error('invalid simulationMethod %s',simulationMethod);
end

end % SimulationGM


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
