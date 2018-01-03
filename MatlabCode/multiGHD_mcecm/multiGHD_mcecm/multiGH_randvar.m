function f=multiGH_randvar(n,mu,lambda,gamma,Sigma,chi,psi)
  
  %
  % [F]=MULTIGH_RANDVAR(N, MU, LAMBDA, GAMMA, SIGMA, CHI, PSI)
  % Draw samples from a Multivariate (multidimensional) Generalized
  % Hyperbolic Distribution (M-GHD). 
  % 
  % N is the number of samples to be drawn.
  %
  % LAMBDA, GAMMA, SIGMA, CHI, and PSI are parameters of a M-GHD. By now
  % you should have got used to these little creatures ;-)
  %
  % F contains the samples drawn
  %
  %
  % NOTE: This code uses the 'randraw' Matlab program. 
  % 'randraw' is available on the Matlab file exchange webpage:
  % http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7309&objectType=File
  %
  %
  % -------------------------------------------------------------------
  % Author : Saket Sathe
  % Email : saket@ee.iitb.ac.in
  % Date : 9th June 2006
  % -------------------------------------------------------------------
  % 

    
  f=[];
  d=max(size(mu));
  norm_dn = mvnrnd(zeros(1,d),Sigma,n);
  w=randraw('gig',[lambda,chi,psi],n);
  
  for i=1:n
    temp = mu + sqrt(w(i))*norm_dn(i,:);
    f = [f;temp];
  end
