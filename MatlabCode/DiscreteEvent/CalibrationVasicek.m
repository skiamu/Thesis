function [param] = CalibrationVasicek(Returns,dt,CalibrationMethod)
%CalibrationVasicek is a function for calibrating the Vasicek model for the
%short-rate via a least-squares method.
%   INPUT:
%      Returns =
%      CalibrationMethod =
%
%   OUTPUT:

switch CalibrationMethod
	case 'LS'
		X = Returns(1:end-1); % regressor
		y = Returns(2:end); % response
		figure
		plot(X,y,'.')
		xlabel('r_{i-1}')
		ylabel('r_i')
		print('/home/andrea/Thesis/Latex/final/Images/LinearRegression','-dpng', '-r900');
		mdl = fitlm(X,y) % 
		Coeff = mdl.Coefficients.Estimate; % [beta, alpha]
		param.a = (1 - Coeff(2)) / dt;
		param.b = Coeff(1) / (param.a * dt);
		param.sigma = sqrt(mdl.MSE / dt);
		param.RSEa = (1 - sqrt(mdl.CoefficientCovariance(2,2))) / dt;
		param.RSEb = sqrt(mdl.CoefficientCovariance(1,1)) / (param.RSEa * dt);
	case 'MLE'
		n = length(Returns);
		Sx = sum(Returns(1:end-1));
		Sy = sum(Returns(2:end));
		Sxx = sum(Returns(1:end-1).^2);
		Sxy = sum(Returns(1:end-1) .* Returns(2:end));
		Syy = sum(Returns(2:end).^2);
		param.b = (Sy * Sxx - Sx * Sxy) / (n * (Sxx - Sxy) - (Sx^2 - Sx * Sy));
		param.a = -log((Sxy - param.b * Sx - param.b * Sy + n * param.b^2) / ...
			(Sxx - 2 * param.b * Sx + n * param.b^2)) / dt;
		alpha = exp(-param.a * dt);
		sigma_tilde_2 = (Syy - 2 * alpha * Sxy + alpha^2 * Sxx - 2 * param.b * ...
			(1 - alpha) * (Sy - alpha * Sx) + n * param.b^2 * (1 -alpha)^2) / n;
		param.sigma = sqrt(sigma_tilde_2 * 2 * param.a / (1 - alpha^2));
	otherwise
		error('invalid CalibrationMethod %s',CalibrationMethod)
end
end % CalibrationVasicek

