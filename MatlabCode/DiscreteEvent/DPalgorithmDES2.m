function [U,J] = DPalgorithmDES2(N,X,p,lambda,r,J_jump,VaR,alpha)
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
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
lb = 0; ub = 1; % upper and lower bound weigth cash
eta = 1e-4 / 3; % integretion interval discretization step
%% optimization
for k = N : -1 : 1
	k % print current iteration
	u0 = 0.5; % initial condition
	dimXk = length(X{k}); % number of single optimizations
	Uk = zeros([dimXk 1]);
	Jk = zeros([dimXk 1]);
	int_domain = (X{k+1}(1):eta:X{k+1}(end))'; % integretion domain with more points
	Jinterp = interp1(X{k+1},J{k+1},int_domain);
	for j = dimXk : -1 : 1
		[Uk(j),Jk(j)] = fmincon(@(u) -objfun(u,X{k}(j),int_domain,Jinterp,p,lambda,r,J_jump),...
			u0,[],[],[],[],lb,ub,@(u)confuneq(u,p,lambda,r,J_jump,VaR,alpha),options);
		u0 = Uk(j);
	end
	U{k} = Uk; J{k} = -Jk;
	% print allocation maps
	if k ~= 1
		idx = find(X{k} <= 1.4 & X{k} >= 0.6);
		figure
		plot(X{k}(idx),U{k}(idx),'b.--')
		title(strcat('k = ',num2str(k-1)))
		legend('cash','stock')
		% 	saveas(gcf,[pwd strcat('/Latex/secondWIP/k',num2str(k-1),model,'.png')]);
	end
end

end


function f = objfun(u,x,int_domain,J,p,lambda,r,J_jump)
%objfun defines the problem's objective function that will be maximized
%   INPUT:
%      u = asset allocation vector
%      x = realization portfolio value [scalar]
%      int_domain = integration domain, k+1 target set discretized
%      J = k+1-th optimal value function
%      param = density parameters
%   OUTPUT:
%      f = objective function
	f = trapz(int_domain, J .* pfDES(int_domain,x,u,J_jump,p,lambda,r));
end % objfun



function [c,ceq] = confuneq(u,p,lambda,r,J_jump,VaR,alpha)
%confuneq states the max problem's constraints
%   INPUT:
%      u = asset allocation vector
%      p = 
%      lambda = 
%      r = 
%      J_jump = 
%      VaR = value at risk
%      alpha = confidence level
%   OUTPUT:
%      c = inequality constraints
%      ceq = equality constraints

Var_w = u^2 * lambda * r^2 / ((lambda - 2*r) * (lambda - r)^2) + ...
	((1-u)*J_jump)^2 * 4 * p * (1-p);
sigmaMax = VaR * sqrt(12) / (sqrt(lambda) * norminv(1-alpha));
c = Var_w - sigmaMax^2;
ceq = [];
end % confuneq






