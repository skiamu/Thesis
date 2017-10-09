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

