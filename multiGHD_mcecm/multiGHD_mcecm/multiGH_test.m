%
% Calibrates a Multivariate Generalized Hyperbolic distribution (M-GHD) for a
% portfolio of assets. This code uses the MCECM based procedure outlined
% in [1].
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

clear all;clc;

orig_x = csvread('your-portfolio-file.csv'); % your-portfolio-file.csv *must* contain the asset prices. One
                                             % asset per column and *must* be comma separated.

T=size(orig_x,1); % row
N=size(orig_x,2); % column

% Compute log returns
for n=1:N
  sing_asset_ret=100*diff(log(orig_x(:,n)));
  new_x(:,n) = sing_asset_ret;
end

clear orig_x; % Get rid of ORIG_X. It is not required anymore.

global debug; % Declare DEBUG and 
debug=0;      % set it.
maxIter=200;  % Maximum number of iterations
tolset = 1e-3; % Tolerance
lambda=-0.5    % Value of LAMBDA at which the Multivariate Generalized
               % Hyperbolic Distribution is estimated.

chi=1          % Initial values of CHI and PSI. A little fiddling and/or
psi=1          % judgement is required to set these values.

% Use MULTIGH_MCECM_CLAM_FIT(). This returns estimates of MU, LAMBDA,
% GAMMA, SIGMA, CHI, PSI for the M-GHD

[mu, lambda, gamma, Sigma, chi,psi]=multiGH_mcecm_clam_fit(new_x,lambda,chi,psi,maxIter,tolset); 


