function [Pout, err] = ob1_calc_power(E, H, mode)
% Description
%     Calculate power in an output mode. Assumes that all the power in the
%     mode is propagating to the right (no reflections).

dims = size(E);

    
    % 
    % Obtain magnitude of mode in both E and H fields.
    %

E_mag = calc_mag(mode.Ey, E);
H_mag = calc_mag(mode.Hz, H);

fields = [E_mag(:), H_mag(:)];

    
    % 
    % Estimate the amplitude, wave-vector and phase of the mode.
    %

phase = unwrap(angle(fields));
x = [[1:dims(1)]-0.5; [1:dims(1)]]'; % Accounts for offset in Yee grid.

A = mean(abs(fields), 1); 
[beta(1), phi(1)] = my_linear_regression(x(:,1), phase(:,1));
[beta(2), phi(2)] = my_linear_regression(x(:,2), phase(:,2));


    %
    % Local optimization to tune fit parameters.
    %

p = mean([A; beta; phi], 2);
fun = @(p) p(1) * exp(i * (p(2) * x + p(3)));
err_fun = @(p) (1/norm(fields, 'fro')) * ...
        norm(fields - fun(p));
options = optimset('MaxFunEvals', 1e3);
p = fminsearch(err_fun, p, options); % Optimize.


    %
    % Calculate output power and error.
    %

amplitude = p(1);
Pout = amplitude^2;

err = err_fun(p);
err_limit = 1e-1;
if (err > err_limit) % If error is somewhat large, tell user.
    warning('Error in fit exceeds threshold (%e > %e).', err, err_limit);
end

% % Plot.
% subplot 111
% plot(real(fields), 'o-');
% hold on
% plot(real(fun(p)), '.-');
% hold off
% pause


function [mag] = calc_mag (ref, field)
% Project y onto x.
my_dot = @(x, y) (dot(y(:), x(:)) / norm(x(:))^2);

for k = 1 : size(field, 1)
    mag(k) = my_dot(ref, field(k,:));
end



function [a, b] = my_linear_regression(x, y)
% Fit data to simple model, y = a * x + b
y_ctr = y - mean(y);
x_ctr = x - mean(x);
a = (x_ctr' * y_ctr) / norm(x_ctr)^2;
b = -(a * mean(x) - mean(y));


