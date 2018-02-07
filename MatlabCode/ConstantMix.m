function [U] = ConstantMix(param,model,VaR,M,N,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% 1) Compute mu and Sigma
switch model
	case 'Gaussian'
		mu = param.mu;
		Sigma = param.S;
	case 'Mixture'
		Sigma = 0;
		mu = 0;
		for i = 1 : length(param)
			mu = mu + param(i).lambda * param(i).mu;
			Sigma = Sigma + param(i).lambda * param(i).S;
			for j = 1 : i-1
				Sigma = Sigma + param(i).lambda * param(j).lambda * ...
					(param(i).mu - param(j).mu) * (param(i).mu - param(j).mu)';
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
		Sigma = wMean * Sigma + wVar * (gamma * gamma');
		mu = param.mu + wMean * gamma;
end

%% 2) Optimization
options = optimoptions(@fmincon,'Algorithm','sqp',...
	'SpecifyConstraintGradient',false);
Nsim = 50;
u0 = rand([Nsim M]);
Aeq = ones([1 M]); beq = 1;
ub = ones([M 1]); lb = zeros([M 1]);
u = zeros([Nsim M]);f = zeros([Nsim 1]);
for i = 1 : Nsim
	[u(i,:),f(i)] = fmincon(@(u)-u' * mu,u0(i,:)',[],[],Aeq,beq,lb,ub,...
		@(u)nonlcon(u,Sigma,VaR,alpha),options);
end
[~,idx] = min(f);
u_star = u(idx,:);
U = cell([N 1]);
for i = 1 : N
	U{i} = u_star;
end

end % ConstantMix


%% 3) non-linear constraint
function [c,ceq] = nonlcon(u,Sigma,VaR,alpha)
sigmaMax = VaR / norminv(1-alpha);
c = sqrt(u' * Sigma * u) - sigmaMax;
ceq = [];
end

