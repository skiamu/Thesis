function [ w ] = SimulationReturns(param,Nsim,M,Nstep,simulationMethod,model,GM )
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
		
	case 'Mixture'
		[ w ] = SimulationGM(GM,param,Nsim,M,Nstep,simulationMethod);
	case 'GH'
		
	otherwise
		error('invalid model %s',model);
end

end

