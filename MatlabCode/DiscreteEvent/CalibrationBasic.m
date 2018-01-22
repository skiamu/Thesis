function [param] = CalibrationBasic(AssetPrice,J_jump,dt)
%CalibrationBasic id a function for calibrating the basic model (parameters p and lambda)
%   INPUT:
%      AssetPrice = risky-asset prices [vector]
%      J_jump = jump size
%   OUTPUT:
%      param = struct of model parameter

%% 1) discrete dymanics
[D,DiscreteS] = DiscretePrice(AssetPrice,J_jump);
figure 
xx = (1:length(D))';
plot(xx,AssetPrice,'b-',xx,DiscreteS,'r-')
legend('Future S&P 500','Discrete Future S&P 500','Location','NorthWest')
xlabel('time')
ylabel('Future S&P 500')
print('/home/andrea/Thesis/Latex/final/Images/DiscreteDynamics',...
	'-dpng', '-r900');


alpha = sum(D == 0);
beta = sum(D == 1);
gamma = sum(D == -1);
param.lambda = -log(alpha/(alpha+beta+gamma))/dt;
param.p = beta / (beta+gamma);
end % CalibrationBasic


