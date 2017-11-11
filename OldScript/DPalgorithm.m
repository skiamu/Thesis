function [ U, J] = DPalgorithm(N,M,X,param,model,VaR,alpha)
%DPalgorithm implements a Dynamic Programming algorithm to solve a
%stochastic reachability problem
%   INPUT:
%      N = number of time steps
%      M = dimension asset allocation (e.g. 3)
%      X = cell array of discretized target sets, to access the i-th
%          element use X{i}
%      param = density parameters (struct or cell array)
%      VaR = yearly value at risk
%      alpha =
%   OUTPUT:
%      U = cell array of asset allocations
%   Remarks:
%      1) when the mixture is used, param is a cell array of the following
%         form {{mu1, Sigma1, lambda1}, ..., {mu_n, Sigma_n, lambda_n}}
%      2) aggiungere VaR nel methodo 1 'Mixture', aggiungere metodo
%      analitico nelle GH

U = cell([N 1]); % cell array, in each position there'll be a matrix
J = cell([N+1 1]); % optimal value function
J{end} = ones([length(X{end}) 1]);  % indicator function of the last target set
options = optimoptions(@fmincon,'Algorithm','active-set','SpecifyConstraintGradient',true,...
	'Display','off','FiniteDifferenceType','central');
% options = optimoptions(@fmincon,'Algorithm','active-set','Display','off','SpecifyConstraintGradient',true);
Aeq = ones([1 M]); beq = 1; % budget constraint
lb = zeros([M 1]); ub = ones([M 1]);
eta = 1e-4; % quantization step for integration
for k = N : -1 : 1
	k % print the current iteration
	u0 = [1; 0.; 0.]; % initial condition
	dimXk = length(X{k});
	Uk = zeros([dimXk M]);
	Jk = zeros([dimXk 1]);
	int_domain = (X{k+1}(1):eta:X{k+1}(end))'; % improve int domain finesse
	Jinterp = interp1(X{k+1},J{k+1},int_domain); % interpolate optimal value function
% 		int_domain = X{k+1};
% 		Jinterp = J{k+1};
	for j = dimXk : -1 : 1
		[Uk(j,:), Jk(j)] = fmincon(@(u)-objfun(u,X{k}(j),int_domain,Jinterp,param,model), ...
			u0,[],[],Aeq,beq,lb,ub,@(u)confuneq(u,param,model,VaR,alpha),options);
 		u0 = Uk(j,:)'; % ititial condition = previuos optial solution
	end
	U{k} = Uk; J{k} = exp(-Jk);
	if k ~= 1
		idx = find(X{k} <= 1.15 & X{k} >= 0.98);
		figure
		area(X{k}(idx),U{k}(idx,:))
		title(strcat('k = ',num2str(k-1)))
		% 	saveas(gcf,[pwd strcat('/Latex/secondWIP/k',num2str(k-1),model,'.png')]);
	end
end
end %  DPalgorithm

function f = objfun(u,x,int_domain,J,param,model)
%objfun defines the problem's objective function that will be maximized
%   INPUT:
%      u = asset allocation vector
%      x = realization portfolio value [scalar]
%      int_domain = integration domain, k+1 target set discretization
%      J = k+1 optimal value function
%      param = density parameters
%   OUTPUT:
%      f = objective function

f = trapz(int_domain, J .* pf(int_domain,x,u,param,model));
f = log(f);
% f = simpsons(J .* pf(int_domain,x,u,param,model),int_domain(1),int_domain(end),[]);
end % objfun

function [c, ceq, D, Deq] = confuneq(u,param,model,VaR,alpha)
%confuneq states the max problem's constraints
%   INPUT:
%      u = asset allocation vector
%      param = density parameters
%      model = 'Gaussian','Mixture','GH'
%      VaR =
%      alpha =
%   OUTPUT:
%      c = inequality constraints
%      ceq = equality constraints
VaRConstraintType = '2';
switch model
	case 'Gaussian'
		S = param.S; mu = param.mu;
		mu_p = -u' * mu; sigma_p = sqrt(u' * S * u);
		ceq = [];
		c = -VaR + (mu_p + norminv(1-alpha) * sigma_p);  % variance constraint
	case 'Mixture'
		n = length(param);
		switch VaRConstraintType
			case '1' % analytic method
				% 1) compute the cdf
				Phi = 0;
				for i = 1 : n
					mu = u' * param{i}{1};
					sigma = sqrt(u' * param{i}{2} * u);
					lambda = param{i}{3};
					Phi = Phi + lambda * normcdf((-VaR - mu) / sigma);
				end
				% 2) setup the constraints
				ceq = [];
				c = Phi - alpha;    % variance constraint
			case '2' % quadratic method
				% 1) compute covariance matrix of asset return w(k+1)
				W = zeros(length(u));
				for i = 1 : n
					W = W + param{i}{3} * param{i}{2};
					for j = 1 : i-1
						W = W + param{i}{3}*param{j}{3}*(param{i}{1}-param{j}{1})...
							*(param{i}{1}-param{j}{1})';
					end
				end
				% 2) setup the constrints
				% By using this formula we suppose gaussianity, mu = 0 and
				% scaling rule, see how it can be improved
				sigma_max = VaR / norminv(1-alpha);
				ceq = [];
				c = u' * W * u - sigma_max^2;
				if nargout > 2
					D = 2 * W * u;
					Deq = [];
				end
		end
		
	case 'GH'
		lambda = param.lambda; Chi = param.Chi; Psi = param.Psi;
		Sigma = param.sigma; gamma = param.gamma;
		% 1) compute Cov[w(k+1)] = E[w(k+1)]*Sigma + Var[w(k+1)]*gamma*gamma'
		wMean = (Chi/Psi)^0.5 * besselk(lambda+1,sqrt(Chi*Psi)) / besselk(lambda,sqrt(Chi*Psi));
		wMoment2 = (Chi/Psi) * besselk(lambda+2,sqrt(Chi*Psi)) / besselk(lambda,sqrt(Chi*Psi));
		wVar = wMoment2 - wMean^2;
		wCov = wMean * Sigma + wVar * (gamma * gamma');
		% 2) set up the constraints
		% By using this formula we suppose gaussianity, mu = 0 and
		% scaling rule, see how it can be improved
		sigma_max = VaR / norminv(1-alpha);
		ceq = [];
		c = u' * wCov * u - sigma_max^2;    % variance constraint
		if nargout > 2
			D = 2 * wCov * u;
			Deq = [];
		end
	otherwise
		error('invalid model %s',model);
end

end % confuneq
















