function f = l2max_mcecm_clam(x,y);

  %
  % F = L2MAX_MCECM_CLAM(X,Y)
  % Function to me maximized with respect to CHI and PSI.
  % Part of the MCECM algorithm with constant LAMBDA [1].
  %
  %
  %
  % References: 
  %       [1] McNeil, A. and Frey, R. and Embrechts, P. (2005) 
  %        Quantitative Risk Management, Princeton University Press.
  %
  % -------------------------------------------------------------------
  % Author : Saket Sathe
  % Email : saket@ee.iitb.ac.in
  % Date : 9th June 2006
  % -------------------------------------------------------------------
  % 

    
  chi=x(1); %  CHI
  psi=x(2); %  PSI
  
  lambda=y(1); % LAMBDA 
  delta=y(2);  % DELTA
  eta=y(3);    % ETA
  n=y(4); 
  xi=y(5);     % XI
  
 f1 = (lambda-1)*xi;
  
 f2 = -0.5*chi*n*delta - 0.5*psi*n*eta - 0.5*n*lambda*log(chi) + 0.5*n*lambda*log(psi);
  
 f3 = - n*log(2*besselk(lambda,sqrt(chi*psi)));
  
 f = - (  (f1) + (f2) + (f3) );

 % f = - ( f2 + f3 ); % Just in case you wish to exlude f1=(lambda-1)*xi
