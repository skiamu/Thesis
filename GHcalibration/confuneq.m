function [A,b,Aeq,beq] = confuneq(GHmodel,d)
%confuneq is a function for setting the constraint for the maximization
%problem
%   INPUT:
%   OUTPUT:

switch GHmodel
	case 'NIG' % lambda = -0.5, alpha_bar > 0
		A = [0 -1]; b = -eps;
		Aeq = [1 0]; beq = -0.5;
	case 'hyp'% lambda = (d+1)/2, alpha_bar > 0
		A = [0 -1]; b = -eps;
		Aeq = [1 0]; beq = (d+1)/2;
	case 'VG'
		% lambda > 0
		A = -1; b = -eps;
		Aeq = []; beq = [];
	case 't'
		% lambda < -1
		Aeq = []; beq = [];
		A = 1; b = -eps -1;
	otherwise % alpha > 0
		Aeq = []; beq = [];
		A = [0 -1]; b = -eps;
end

end % objfun
