function [U,J] = ODAAalgorithmDES(N,X,J_jump,param,model)
%DPalgorithm implements a Dynamic Programming algorithm to solve a
%stochastic reachability problem
%   INPUT:
%      N = number of events
%      M = dimension asset allocation (e.g. )
%      X = cell array of discretized target sets, to access the i-th
%          element use X{i}
%      J_jump = jump size
%      param = struct of parameters
%   OUTPUT:
%      U = cell array of asset allocations
%      J = cell array of optial value functions

%% initialization
U = cell([N 1]); % asset allocation cell array
J = cell([N+1 1]); % optimal value function cell  array
J{end} = ones([length(X{end}) 1]); % indicator function target set X_N
options = optimoptions(@fmincon,'Algorithm','active-set','Display','off');
lb = -1; ub = 1; % short positions on the risky asset are allowed
eta = 1e-4/3; % integretion interval discretization step (1e-4/2 for ext1)
if strcmp(model,'ext1') % risk constraint ext1
	global VarTau ExpValue_tau;
	mu = param.mu; sigma=param.sigma;mu_tilde = mu-0.5*sigma^2;
	a = @(n) 0.5 * mu_tilde^2 / sigma^2 + sigma^2 * n.^2 * pi^2 / (8 * J_jump^2);
	epsilon = 1e-8;
	N1 = ceil(((8*J_jump^2/(sigma^2*pi^2))^3/(epsilon*4))^(1/4));
	NN1 = 1:N1;
	N2 = ceil(((8*J_jump^2/(sigma^2*pi^2))^2/(epsilon*2))^(1/2));
	NN2 = 1:N2;
	ExpValue_tau2 = 2*cosh(mu_tilde*J_jump/sigma^2)*sigma^2*pi/(4*J_jump^2)*...
		sum(2*NN1.*(-1).^(NN1+1)./a(NN1).^3 .*sin(NN1*pi/2));
	ExpValue_tau = 2*cosh(mu_tilde*J_jump/sigma^2)*sigma^2*pi/(4*J_jump^2)*...
		sum(NN2.*(-1).^(NN2+1)./a(NN2).^2 .*sin(NN2*pi/2));
	VarTau = ExpValue_tau2 - ExpValue_tau^2;
end
%% optimization
for k = N : -1 : 1
	k % print current iteration
	dimXk = length(X{k}); % number of single optimizations
	Uk = zeros([dimXk 1]);
	Jk = zeros([dimXk 1]);
	int_domain = (X{k+1}(1):eta:X{k+1}(end))'; % integretion domain with more points
	Jinterp = interp1(X{k+1},J{k+1},int_domain);
	if k ~= 1
		u0 = -1; % initial condition
	else
		u0 = 1;
	end
	for j = dimXk:-1:1
		j % print current iteration
		[Uk(j),Jk(j)] = fmincon(@(u) -objfun(u,X{k}(j),int_domain,Jinterp,param,J_jump,model),...
			u0,[],[],[],[],lb,ub,[],options);
		u0 = Uk(j)
	end
	U{k} = Uk; J{k} = -Jk;
	% print allocation maps
	if k ~= 1
		idx = find(X{k} <= 3 & X{k} >= 0.1);
		figure
		area(X{k}(idx),U{k}(idx))
		title(strcat('k = ',num2str(k-1)))
		% 	saveas(gcf,[pwd strcat('/Latex/secondWIP/k',num2str(k-1),model,'.png')]);
	end
end

end


function f = objfun(u,x,int_domain,J,param,J_jump,model)
%objfun defines the problem's objective function that will be maximized
%   INPUT:
%      u = asset allocation vector
%      x = realization portfolio value [scalar]
%      int_domain = integration domain, k+1 target set discretized
%      J = k+1-th optimal value function
%      J_jump = jump size
%      param = struct of parameters
%   OUTPUT:
%      f = objective function
switch model
	case 'basic'
		f = trapz(int_domain, J .* pfDES(int_domain,x,u,J_jump,param));
	case 'ext1'
		f = trapz(int_domain, J .* pfDESext1(int_domain,x,u,J_jump,param));
	case 'ext2'
		f = trapz(int_domain', J' .* pfDESext2(int_domain',x,u,J_jump,param));
end
end % objfun


function [c, ceq] = constr(u,param,J_jump,model)
alpha = 0.01;
VaR_m = 0.07;
if strcmp(model,'basic')
	r = param.r; lambda = param.lambda; p = param.p;
	sigma_max = sqrt(12) * VaR_m  / (sqrt(lambda) * norminv(1-alpha));
	c = lambda * r^2 / ((lambda-2*r)*(lambda-r)^2) + (u*J_jump)^2 *...
		4*p*(1-p) - sigma_max^2;
else % ext1 model
	global VarTau ExpValue_tau;
	r = param.r;
	p = param.p;
	sigma_max = VaR_m * sqrt(12*ExpValue_tau) / norminv(1-alpha);
	c = r^2*VarTau + u^2*4*J_jump^2*p*(1-p) - sigma_max^2;
end
ceq = [];
end



