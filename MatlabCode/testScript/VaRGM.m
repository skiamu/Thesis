function y = VaRGM(param,u,alpha)
y0 = 0.02;
y = zeros([1 length(alpha)]);
options = optimset('Display','off');
for i = 1 : length(alpha);
	y(i) = fsolve(@(z) Phi(param,u,z,alpha(i)),y0,options);
	y0 = y(i);
end
end % ValueAtRisk

function y = Phi(param,u,z,alpha)
Phi = 0;
for i = 1 : length(param)
	mu = -u' * param(i).mu;
	sigma = sqrt(u' * param(i).S * u);
	Phi = Phi + param(i).lambda * normcdf(-(z - mu) / sigma);
end
y = Phi - alpha;
end