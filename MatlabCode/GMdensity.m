
function f = GMdensity(z,param,k)
%GMdensity is the density of a random vector that follows a gaussian
%mixture distribution
%   INPUT: 
%      z = point where to compute the density, if z is a matrix the density
%          is computed at each row [array or matrix]
%      param = struct array or parameters
%      k = number of gaussian components
%   OUTPUT:
%      f = density value at z
f = 0;
for i = 1 : k
	f = f + param(i).lambda * mvnpdf(z,param(i).mu',param(i).S);
end
end % GMdensity

