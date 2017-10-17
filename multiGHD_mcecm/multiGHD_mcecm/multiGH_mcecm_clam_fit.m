function [mu, lambda, gamma, Sigma, chi, psi]=multiGH_mcecm_clam_fit(X,lambda,chi,psi,maxIter,tolset)

  %
  % [MU, LAMBDA, GAMMA, Sigma, CHI, PSI]=MULTIGH_MCECM_CLAM_FIT(X, LAMBDA, MAXITER, TOLSET)
  % Calibration of a Multivariate (multidimensional) Generalized
  % Hyperbolic Distribution (M-GHD) using the Multi-cycle Expectation Conditional
  % Maximization (MCECM) algorithm for constant value of LAMBDA [1].
  %
  % X is a matrix of random variables for which a M-GHD is to be
  % calibrated. X=[X_1 X_2 ...X_N], where each column X_N is of length T.
  % Each row of X contains samples form a N-dimensional probability distribution.
  %
  % LAMBDA is the constant value at which the M-GHD is estimated.
  % 
  % CHI and PSI are their initial values. A little fiddling and/or
  % judgement is required to set these values.
  %
  % MAXITER is the maximum number of iterations the algorithm should iterate.
  %
  % TOLSET is the tolerence for detecting convergence of algorithm. 
  % The algorithm is declared converged and quit if the maximum absolute
  % value of the percent change in estimated parameters is < TOLSET
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

  
    global debug; % DEBUG is a global variable. It can be toggled for
                  % debug messages. DEBUGMSG() is used to display the
                  % debug messages
    
    T=size(X,1); % row
    N=size(X,2); % column
  
    options = optimset('Display','notify','TolFun',1e-5,'MaxFunEvals',5000, ...
                       'MaxIter',2000);
    
    % Initial Values
    
    mu = mean(X);
    gamma = zeros(1,N);
    Sigma = cov(X);
    S = Sigma;
    chi = chi; % Redundant ;-)
    psi = psi;
    
    %Tolerance
    
    tol = tolset;
    
    
    if debug==1
      disp('Initial Values:');
      debugmsg('lambda:',lambda,0);
      debugmsg('mu:',mu,0);
      debugmsg('gamma:',gamma,0);
      debugmsg('Sigma:',Sigma,0);
      debugmsg('chi:',chi,0);
      debugmsg('psi:',psi,1);
    end
    
    chis=[];
    psis=[];
    
    for iter=1:maxIter % Main iteration
      disp('----------------------------------------------------');
      disp(sprintf('Iteration: %d',iter));
      disp('----------------------------------------------------');
      
      delta_i=zeros(1,T);
      eta_i=zeros(1,T);
      xi_i=zeros(1,T);
      
      for i=1:T
        rho = (X(i,:)-mu)*inv(Sigma)*(X(i,:)-mu)';
        f1=((rho+chi)/(psi+gamma*inv(Sigma)*gamma')) ;
        f2=sqrt((rho+chi)*(psi+gamma*inv(Sigma)*gamma'));
        
        delta_i(i)=((f1)^(-0.5))*((besselk(lambda-(N/2)-1,f2))/(besselk(lambda-(N/2),f2)));
        eta_i(i)=((f1)^(0.5))*((besselk(lambda-(N/2)+1,f2))/(besselk(lambda-(N/2),f2)));
        xi_i(i)=0.5*log(f1);
      end
      
      delta=(1/T)*sum(delta_i); % delta[k]
      eta=(1/T)*sum(eta_i);     % eta[k]
      xi=(1/T)*sum(xi_i);       % xi[k]
      
      debugmsg('delta value:',delta,1);
      debugmsg('eta value:',eta,1);
      debugmsg('xi value:',xi,1);
      
      old_params=[eta delta xi gamma mu chi psi];
      for i=1:size(Sigma)
        old_params=[old_params Sigma(i,:)];
      end
      debugmsg('old_params',old_params,1);
      
      gamma = zeros(1,N);
      for i=1:T
        gamma=gamma +delta_i(i)*(mu-X(i,:));
      end
      gamma=((1/T)*gamma)/(delta*eta -1); % gamma[k+1]
      
      debugmsg('gamma value:',gamma,1);
      
      mu = zeros(1,N);
      for i=1:T
        mu=mu +delta_i(i)*X(i,:);
      end
      mu = ((1/T)*mu-gamma)/(delta); %mu[k+1]
      
      debugmsg('mu value:',mu,1);
      
      Psi=zeros(N,N);
      for i=1:T
        Psi = Psi + delta_i(i)*(X(i,:)-mu)'*(X(i,:)-mu);
      end
      
      Psi = (1/T)*Psi-eta*gamma'*gamma; 
      Sigma = ((det(S)^(1/N))/(det(Psi)^(1/N)))*Psi; % Sigma[k+1]
      
      debugmsg('Psi value:',Psi,1);
      debugmsg('Sigma value:',Sigma,1);
      
      delta_i=zeros(1,T);
      eta_i=zeros(1,T);
      xi_i=zeros(1,T);
      for i=1:T
        rho = (X(i,:)-mu)*inv(Sigma)*(X(i,:)-mu)';
        f1=((rho+chi)/(psi+gamma*inv(Sigma)*gamma')) ;
        f2=sqrt((rho+chi)*(psi+gamma*inv(Sigma)*gamma'));
        
        delta_i(i)=((f1)^(-0.5))*((besselk(lambda-(N/2)-1,f2))/(besselk(lambda-(N/2),f2)));
        eta_i(i)=((f1)^(0.5))*((besselk(lambda-(N/2)+1,f2))/(besselk(lambda-(N/2),f2)));
        xi_i(i)=0.5*log(f1);
      end
      
      delta=(1/T)*sum(delta_i); % delta[k,2]
      eta=(1/T)*sum(eta_i);     % eta[k,2]
      xi=(1/T)*sum(xi_i);       % xi[k]
      
      %
      %  [estimates fval exflag output] = fminsearch(@l2max,[chi psi],options,[lambda ...
      %                      delta eta T xi]);
      
      a = [-1 0; 0 -1];
      b= [-0.01;-0.01];
      [estimates fval exflag output] = fmincon(@l2max_mcecm_clam,[chi, psi],a,b,[],[],[],[],[],options,[lambda ...
                          delta eta T xi]);  % Call L2MAX_MCECM_CLAM() for constrained
                                             % minimization of the likelihood function to
                                             % estimate CHI and PSI.
      
      chi=estimates(1) % Pick-up fresh estimate of CHI
      psi=estimates(2) % Pick-up fresh estimate of PSI
      
      chis = [chis chi]; % Push CHI in an array. Used for debugging.
      psis = [psis psi]; % Push PSI in an array. Used for debugging.
      
      new_params=[eta delta xi gamma mu chi psi];
      
      for i=1:size(Sigma)
        new_params=[new_params Sigma(i,:)];
      end
      
      paramdiff = new_params-old_params;
      perdiff = paramdiff./old_params;
      
      debugmsg('difference in parameters:',paramdiff,1);
      debugmsg('maximum absolute value of the percentage change:',max(abs(perdiff)),1);
      
      if max(abs(perdiff))<tol
        disp('I am done. Finally! ;-)');
        break
      end
      %  pause % Used for debugging
    end
