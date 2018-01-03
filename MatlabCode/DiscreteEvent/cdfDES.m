function [ F ] = cdfDES( z,x,u,J,p,lambda,r )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
F1 = zeros(size(z));
F2 = F1;
csi = x * J * u;
idx1 = x - csi < z & z < x + csi;
idx2 = z >= x + csi;
F1(idx1) =  (1 - p) * (1 - ((z(idx1) + csi) / x).^(-(lambda + 1) / r));
F2(idx2) = (1-p) * (1 - ((z(idx2) + csi) / x).^(-(lambda + 1) / r)) + ...
	p * (1 - ((z(idx2) - csi) / x).^(-(lambda + 1) / r));
F = F1 + F2;
end

