function Q = quadvgk(fv, Subs, NF)
%QUADVGK: Vectorised adaptive G7,K15 Gaussian quadrature on vector of integrals
%
%	This function calculates the integration of a vector of functions via adaptive G7,K15
%	Gaussian quadrature. The procedure follows that of Shampine [1]. This function is
%	similar to quadv, but uses the quadgk algorithm.
%
%Usage:
%	Q = quadvgk(fv, Subs, NF)
%
%	Q		= Returned vector of numerical approximations to the integrals
%	fv		= Function handle which returns a matrix of values (see below)
%	Subs	= Matrix of intervals (see below)
%	NF		= Number of functions to be calculated (see below)
%
%Example:
%	Y = @(x,n) 1./((1:n)+x);
%	Qv = quadvgk(@(x) Y(x,10), [0;1], 10);
%
%	Y is a vector of functions, where each row represents a different function and each
%	column corresponds to a different point where the function is evaluated. The function
%	needs to return a matrix of size (NF, NX) where NF is the number of functional
%	integrals to evaluate, and NX is the number of data points to evaluate. In the above
%	example, the quadvgk calculates:
%	[int(1./(1+x), 0..1); int(1./(2+x), 0..1); int(1./(3+x), 0..1); ...; int(1./(10+x), 0..1)]
%
%	Subs is the initial set of subintervals to evaluate the integrand over. Should be in
%	the form [a1 a2 a3 ... aN; b1 b2 b3 ... bN], where "an" and "bn" are the limits of
%	each subinterval. In the above example, Subs = [a; b]. If the function contains sharp
%	features at known positions, the boundaries of these features should be added to Subs.
%	Note that in order to integrate over a continuous span of subintervals, Subs = [A a1
%	a2 a3 ... aN; a1 a2 a3 ... aN B] where "A" and "B" are the limits of the whole
%	integral, and a1..aN are points within this range.
%
%   Based on "quadva" by Lawrence F. Shampine.
%   Ref: L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab",
%   Journal of Computational and Applied Mathematics, to appear.

%	Positions and weightings for Gaussian Quadrature
	%	15-point Kronrod positions & weightings
	XK = [-0.991455371120813; -0.949107912342759; -0.864864423359769; -0.741531185599394; ...
		-0.586087235467691; -0.405845151377397; -0.207784955007898; 0; ...
		0.207784955007898; 0.405845151377397; 0.586087235467691; ...
		0.741531185599394; 0.864864423359769; 0.949107912342759; 0.991455371120813];
	WK = [0.022935322010529; 0.063092092629979; 0.104790010322250; 0.140653259715525; ...
		0.169004726639267; 0.190350578064785; 0.204432940075298; 0.209482141084728; ...
		0.204432940075298; 0.190350578064785; 0.169004726639267; ...
		0.140653259715525; 0.104790010322250; 0.063092092629979; 0.022935322010529];

	%	7-point Gaussian weightings
	WG = [0.129484966168870; 0.279705391489277; 0.381830050505119; 0.417959183673469; ...
		0.381830050505119; 0.279705391489277; 0.129484966168870];
	NK = length(WK);
	G = (2:2:NK)';		%	7-point Gaussian poisitions (subset of the Kronrod points)

	Q = zeros(NF, 1);
	while (~isempty(Subs))
		GetSubs;
		M = (Subs(2,:)-Subs(1,:))/2;
		C = (Subs(2,:)+Subs(1,:))/2;
		NM = length(M);
		x = reshape(XK*M + ones(NK,1)*C, 1, []);
		FV = fv(x);
		Q1 = zeros(NF, NM);
		Q2 = zeros(NF, NM);
		for n=1:NF
			F = reshape(FV(n,:), NK, []);
			Q1(n,:) = M.*sum((WK*ones(1,NM)).*F);
			Q2(n,:) = M.*sum((WG*ones(1,NM)).*F(G,:));
		end
		ind = find((max(abs((Q1-Q2)./Q1), [], 1)<=1e-6)|((Subs(2,:)-Subs(1,:))<=eps));
		Q = Q + sum(Q1(:, ind), 2);
		Subs(:, ind) = [];
	end

	function GetSubs
		M = (Subs(2,:)-Subs(1,:))/2;
		C = (Subs(2,:)+Subs(1,:))/2;
		I = XK*M + ones(NK,1)*C;
		A = [Subs(1,:); I];
		B = [I; Subs(2,:)];
		Subs = [reshape(A, 1, []); reshape(B, 1, [])];
	end
end
	