function f = objfun(x,delta,eta,csi,GHmodel)
%objfun is the objective function Q2 for the maximization problem
%   INPUT:
%      x = vector of parameters, [lambda,Chi,Psi]
%      delta =
%      eta =
%      csi =
%   OUTPUT:
%      f =
% REMARKS: different objective function for 't' and 'VG'
[N, ~] = size(delta);
switch GHmodel
	case 't'
		nu = -2 * x; Chi = nu - 2;
		f = -x * N * log(Chi / 2) - N * log(gamma(-x)) + (x-1) * sum(csi) - ...
			0.5 * Chi * sum(delta);
	case 'VG'
		Psi = 2 * x;
		f = x * N * log(Psi / 2) - N * log(gamma(x)) + (x - 1) * sum(csi) - ...
			0.5 * Psi * sum(eta);
	otherwise % NIG, hyperbolic, generalized hyperbolic
		% x(1) = lambda, x(2) = alpha_bar
		Psi = x(2) * besselk(x(1) + 1,x(2)) / besselk(x(1),x(2));
		Chi = x(2)^2 / Psi;
		f = (x(1)-1) * sum(csi) - 0.5 * Chi * sum(delta) - 0.5 * Psi * sum(eta) -...
			0.5 * N * x(1) * log(Chi) + 0.5 * N * x(1) * log(Psi) - ...
			N * log(2 * besselk(x(1),x(2)));
end

f = - f; % maximization problem
end % objfun