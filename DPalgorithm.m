function [U,J] = DPalgorithm(N,M,X,param,model,VaR,alpha)
%DPalgorithm implements a Dynamic Programming algorithm to solve a
%stochastic reachability problem
%   INPUT:
%      N = number of time steps
%      M = dimension asset allocation (e.g. 3)
%      X = cell array of discretized target sets, to access the i-th
%          element use X{i}
%      param = density parameters [struct array]
%      VaR = value at risk (horizon according to the returns)
%      alpha = confidence level
%   OUTPUT:
%      U = cell array of asset allocations
%      J = cell array of optial value functions

%% initialization
U = cell([N 1]); % asset allocation cell array
J = cell([N+1 1]); % optimal value function cell  array
J{end} = ones([length(X{end}) 1]); % indicator function target set X_N
options = optimoptions(@fmincon,'Algorithm','active-set','SpecifyConstraintGradient',true,...
	'FiniteDifferenceType','central','Display','off');
Aeq = ones([1 M]); beq = 1; % equality constraint
lb = zeros([M 1]); ub = ones([M 1]); % upper and lower bound
eta = 1e-4 / 2; % integretion interval discretization step
% compute covariance matrix
%% optimization
for k = N : -1 : 1
	k % print current iteration
	u0 = [1; 0; 0]; % initial condition
	dimXk = length(X{k}); % number of single optimizations
	Uk = zeros([dimXk M]);
	Jk = zeros([dimXk 1]);
	int_domain = (X{k+1}(1):eta:X{k+1}(end))'; % integretion domain with more points
	Jinterp = interp1(X{k+1},J{k+1},int_domain);
	Risk = zeros([dimXk 1]);
	for j = dimXk : -1 : 1
		[Uk(j,:),Jk(j)] = fmincon(@(u) -objfun(u,X{k}(j),int_domain,Jinterp,param,model,k),...
			u0,[],[],Aeq,beq,lb,ub,@(u)confuneq(u,param,model,VaR,alpha),options);
		[Risk(j),~] = confuneq(Uk(j,:)',param,model,VaR,alpha);
		u0 = Uk(j,:)';
	end
	U{k} = Uk; J{k} = -Jk;
	% print allocation maps
	if k ~= 1
		idx = find(X{k} <= 1.4 & X{k} >= 0.6);
		figure
		area(X{k}(idx),U{k}(idx,:))
		title(strcat('k = ',num2str(k-1)))
		% 	saveas(gcf,[pwd strcat('/Latex/secondWIP/k',num2str(k-1),model,'.png')]);
	end
end

end


function f = objfun(u,x,int_domain,J,param,model,k)
%objfun defines the problem's objective function that will be maximized
%   INPUT:
%      u = asset allocation vector
%      x = realization portfolio value [scalar]
%      int_domain = integration domain, k+1 target set discretized
%      J = k+1-th optimal value function
%      param = density parameters
%   OUTPUT:
%      f = objective function
% if k == 52
% 	mu1 = x * (1 + u' * param(1).mu);
% 	mu2 = x * (1 + u' * param(2).mu);
% 	sigma1 = sqrt(x^2 * u' * param(1).S * u);
% 	sigma2 = sqrt(x^2 * u' * param(2).S * u);
% 	f = param(1).lambda * (1 - normcdf((int_domain(1)-mu1) / sigma1)) + ...
% 		param(2).lambda * (1 - normcdf((int_domain(1)-mu2) / sigma2));	
% else
	f = trapz(int_domain, J .* pf(int_domain,x,u,param,model));
end % objfun



function [c,ceq, D, Deq] = confuneq(u,param,model,VaR,alpha)
%confuneq states the max problem's constraints
%   INPUT:
%      u = asset allocation vector
%      param = density parameters
%      model = 'Gaussian','Mixture','GH'
%      VaR = value at risk
%      alpha = confidence level
%   OUTPUT:
%      c = inequality constraints
%      ceq = equality constraints

switch model
	case 'Gaussian'
		c = -u' * param.mu + norminv(1-alpha) * sqrt(u' * param.S * u) - VaR;
		ceq = [];
	case 'Mixture'
		n = length(param);
		method = '';
		if strcmp(method,'analitic'); % analitic method
			Phi = 0;
			for i = 1 : n
				mu = -u' * param(i).mu;
				sigma = sqrt(u' * param(i).S * u);
				Phi = Phi + param(i).lambda * normcdf(-(VaR - mu) / sigma);
			end
			c = Phi - alpha;
			ceq = [];
		else % quadratic method
			W = 0;
			for i = 1 : length(param)
				W = W + param(i).lambda * param(i).S;
				for j = 1 : i-1
					W = W + param(i).lambda * param(j).lambda * ...
						(param(i).mu - param(j).mu) * (param(i).mu - param(j).mu)';
				end
			end
			sigmaMax = VaR / norminv(1-alpha);
			c = sqrt(u' * W * u) - sigmaMax;
			ceq = [];
			if nargout > 2
				D =  (W * u) / sqrt(u' * W * u);
				Deq = [];
			end
		end
		
	case {'GH','H','NIG','t'}
		if strcmp(model,'t')
			% 0) extract parameters
			nu = param.nu; lambda = param.lambda;
			Sigma = param.sigma; gamma = param.gamma;
			% 1) compute Cov[w(k+1)] = E[w(k+1)]*Sigma + Var[w(k+1)]*gamma*gamma'
			Chi = nu  - 2;
			wMean = 0.5 * Chi / (-lambda - 1);
			wVar = 0.25 * Chi^2 / ((-lambda-1)^2 * (-lambda-2));
		else
			% 0) extract parameters
			lambda = param.lambda; Chi = param.Chi; Psi = param.Psi;
			Sigma = param.sigma; gamma = param.gamma;
			% 1) compute Cov[w(k+1)] = E[w(k+1)]*Sigma + Var[w(k+1)]*gamma*gamma'
			wMean = (Chi/Psi)^0.5 * besselk(lambda+1,sqrt(Chi*Psi)) / besselk(lambda,sqrt(Chi*Psi));
			wMoment2 = (Chi/Psi) * besselk(lambda+2,sqrt(Chi*Psi)) / besselk(lambda,sqrt(Chi*Psi));
			wVar = wMoment2 - wMean^2;
		end
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
		error('Invalid model : %s', model);
end

end % confuneq






