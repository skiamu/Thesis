function [param] = HorizonCorrection(freq,param,model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch freq
	case 'm'
		t = 4;
	case 'q'
		t = 4 * 3;
	case 's'
		t = 4 * 6;
	case 'y'
		t = 52;
	otherwise
		return
end

% we approximate linear returns with log returns
switch model
	case 'Gaussian'
		param.mu = param.mu * t;
		param.S = param.S * t;
	case 'Mixture'
		paramNew(t).mu = 0;
		for j = 0 : t
			paramNew(j+1).mu = param(1).mu * (t-j) + param(2).mu * j;
			paramNew(j+1).S = param(1).S * (t-j) + param(2).S * j;
			paramNew(j+1).lambda = nchoosek(t,j) * param(1).lambda^(t-j) * ...
				param(2).lambda^j;
		end
		param = paramNew;
	case {'GH','NIG','t'}
		param.Chi = param.Chi * t^2 ;
		param.mu = param.mu * t;
end
	
end

