function [a_fit, r_fit, resnorm] = fitmodel_lm(y,TE, params_initial)
% levenberg-marquardt

fun = @(x,xdata)x(1)*exp(-x(2)* xdata);
% Initialize the parameters [a, r]
% p(1) = a, p(2) = r

num_pixels = size(y, 2);
a_fit = zeros(num_pixels, 1);
r_fit = zeros(num_pixels, 1);
resnorm = zeros(num_pixels, 1);

% Create options structure for Levenberg-Marquardt algorithm
options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');

% Fit the model to the data for each pixel
for i = 1:num_pixels
    [pfit, resnorm(i), ~, ~] = lsqcurvefit(fun, params_initial, TE, y(:, i), [], [], options);
    % if exitflag <= 0
    %     a_fit(i) = NaN;  % Indicate non-convergence for this pixel
    %     r_fit(i) = NaN;  
    % else
    %     a_fit(i) = pfit(1);
    %     r_fit(i) = pfit(2);
    % end
    a_fit(i) = pfit(1);
    r_fit(i) = pfit(2);
end

disp(' ------- T2* Estimation using LM Completed ------- ')
end