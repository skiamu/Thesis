% contour plot for different multivariate pdf
clc;close all; clear variables;
addpath(genpath(pwd))

freq = 'd';
M = 3;
[Returns,~,stocks] = getReturns( freq, M );
[param,CalibrationData] = GMcalibrationEM( Returns(:,2:3),2);
%% 1) Gaussian




%% 2) Gaussian Mixture
% 
% Mu = [param(1).mu';param(2).mu'];
% Sigma = cat(3,param(1).S,param(2).S);
% % P = [param(1).lambda,param(2).lambda];
Mu = [0.1 0.8; -.1 -.2];
Sigma = cat(3,[.4 0.3;0.3 .5],[1 -0.9;-0.9 1]);
P = [0.5 0.5];
gm = gmdistribution(Mu,Sigma,P);

gmPDF = @(x,y)pdf(gm,[x y]);
figure;
fcontour(gmPDF,[-2 1.5 -1.5  2],'Fill','on');
hold on
title('GMM - PDF Contours');
saveas(gcf,'/home/andrea/Thesis/Latex/final/Images/GMdensity.png');
figure;
fsurf(gmPDF,[-2 1.5 -1.5  2])
title('PDF of the GMM');

