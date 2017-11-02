find(X{52}== 0.9)
u1 = Uk(2201,:)';
u2 = [0; u1(1); u1(3)];
XN = X{53};
x = 1.1;
figure
plot(XN,pf(XN,x,u1,param,model),'b',XN,pf(XN,x,u2,param,model),'r')
legend('cash','bond')




mu1 = param{1}{1}; sigma1 = param{1}{2}; lambda1 = param{1}{3};
mu2 = param{2}{1}; sigma2 = param{2}{2}; lambda2 = param{2}{3};
lambda = [lambda1 lambda2];
mu = [x * (1 + u1' * mu1),  x * (1 + u1' * mu2)];
sigma = sqrt([x^2 * u1' * sigma1 * u1, x^2 * u1' * sigma2 * u1]);

[a,b,c,d] = ComputeMomentsGM(mu,sigma,[lambda1 lambda2])
W = zeros(length(u1));
for i = 1 : 2
	i
	W = W + param{i}{3} * param{i}{2};
	for j = 1 : i-1
		j
		W = W + param{i}{3}*param{j}{3}*(param{i}{1}-param{j}{1})...
			*(param{i}{1}-param{j}{1})';
	end
end
sqrt(x^2 * u1' * W * u1)

density = @(z,u,x) lambda1 * normpdf(z,x .* (1 + u' * mu1),sqrt(x.^2 * (u'*sigma1*u))) + ...
	lambda2 * normpdf(z,x .* (1 + u' * mu2),sqrt(x.^2 * (u'*sigma2*u)));

u = u1;
mu = [x * (1 + u' * mu1),  x * (1 + u' * mu2)];
sigma = sqrt([x^2 * u' * sigma1 * u, x^2 * u' * sigma2 * u]);

[muC,sigmaC,skC,kurtC] = ComputeMomentsGM(mu,sigma,lambda)
n = length(param); % param is as long as the mixture addendum
f = 0;
u = u1;
for i = 1 : n
	mu = x * (1 + u' * param{i}{1}); % mean f(x,u,w(k+1))
	sigma = sqrt(x^2 * u' * param{i}{2} * u); % std f(x,u,w(k+1))
	lambda = param{i}{3};
	f = f + lambda * normpdf(XN,mu,sigma);
%    muX = muX + lambda * mu;
end
AreaB = trapz(XN,f)
u = u2;
mu = [x * (1 + u' * mu1),  x * (1 + u' * mu2)];
sigma = sqrt([x^2 * u' * sigma1 * u, x^2 * u' * sigma2 * u]);
[muB,sigmaB,skB,kurtB] = ComputeMomentsGM(mu,sigma,lambda)
f = 0;
u = u2;
for i = 1 : n
	mu = x * (1 + u' * param{i}{1}); % mean f(x,u,w(k+1))
	sigma = sqrt(x^2 * u' * param{i}{2} * u); % std f(x,u,w(k+1))
	lambda = param{i}{3};
	f = f + lambda * normpdf(XN,mu,sigma);
%    muX = muX + lambda * mu;
end
AreaC = trapz(XN,f)


figure
xx = 0.9:0.0001:1.06;
plot(xx,density(xx,u1,x),'b',xx,density(xx,u2,x),'r')
legend('cash','bond')

figure
xx = 0.9:0.0001:1.06;
plot(xx,density(x,u1,xx),'b',xx,density(x,u2,xx),'r')
legend('bond','cash')



ComputeMomentsGM(param,u1,x)

u = u2;
mu1 = param{1}{1}; sigma1 = param{1}{2};
mu2 = param{2}{1}; sigma2 = param{2}{2};
muX = 0;
varX = 0;
for i = 1 : n
	mu = x * (1 + u' * param{i}{1}) % mean f(x,u,w(k+1))
	sigma = sqrt(x^2 * u' * param{i}{2} * u) % std f(x,u,w(k+1))
	lambda = param{i}{3};
% 	f = f + lambda * normpdf(z,mu,sigma);
   muX = muX + lambda * mu;
end