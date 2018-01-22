function [U,J] = DPalgorithmDES(N,X,J_jump,param,model)
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
options = optimoptions(@fmincon,'Algorithm','interior-point');
lb = -1; ub = 1; % short positions on the risky asset are allowed
eta = 1e-4/2; % integretion interval discretization step (1e-4/2 for ext1)
%% optimization
for k = N : -1 : 1
	k % print current iteration
	u0 = -1; % initial condition
	dimXk = length(X{k}); % number of single optimizations
	Uk = zeros([dimXk 1]);
	Jk = zeros([dimXk 1]);
	int_domain = (X{k+1}(1):eta:X{k+1}(end))'; % integretion domain with more points
	Jinterp = interp1(X{k+1},J{k+1},int_domain);
	% 	Risk = zeros([dimXk 1]);
	for j = dimXk:-1:1
		j % print current iteration
		[Uk(j),Jk(j)] = fmincon(@(u) -objfun(u,X{k}(j),int_domain,Jinterp,param,J_jump,model),...
			u0,[],[],[],[],lb,ub,[],options);
		u0 = Uk(j);
		% 		[Risk(j),~] = confuneq(u0,p,lambda,r,J_jump,VaR,alpha);
	end
	U{k} = Uk; J{k} = -Jk;
	% print allocation maps
	if k ~= 1
		idx = find(X{k} <= 1.9 & X{k} >= 0.6);
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






