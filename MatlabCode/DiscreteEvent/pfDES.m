function [ f ] = pfDES(z,x,u,J_jump,param)
%pfDES computes the density function of the random variable x(k+1)
%(portdolio value at event number k+1)
%   INPUT:
%      z = indipendent variable
%      x = portfolio value last event
%      u = cash weigth
%      J_jump = jump treshold
%      param = struct of model parameters
%   OUTPUT:
%      f = density computed in z
f = zeros(size(z));

csi = x * J_jump * u;

idx1 = z >= x + csi;

idx2 = z >= x - csi;

p = param.p; lambda = param.lambda; r = param.r;

f(idx1) = p * ((z(idx1) - csi) / x).^(-(lambda + r) / r);

f(idx2) =  f(idx2) + (1-p) * ((z(idx2) + csi) / x).^(-(lambda + r) / r);

f = lambda / (r * x) * f;

end % pfDES




