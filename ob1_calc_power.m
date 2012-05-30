function [Pout, err, filtered_percent] = ob1_calc_power(E, H, mode)
% Description
%     Calculate power in an output mode. Assumes that all the power in the
%     mode is propagating to the right (no reflections).

dims = size(E);

    
    % 
    % Obtain magnitude of mode in both E and H fields.
    %

Pout_0 = mean(sum(real(conj(E) .* H), 2)); % Calculate power in output mode.

E = my_project(mode.Ey, E); % Eliminate other modes from output field.
H = my_project(mode.Hz, H);

Pout = sum(real(conj(E) .* H), 2); % Calculate power in output mode.
err = std(Pout); % Calculate error.
Pout = mean(Pout);
filtered_percent = Pout / Pout_0;

function [f] = my_project(ref, field)
% Project y onto x.

my_dot = @(x, y) dot(y(:), x(:)/norm(x(:))) * x(:) / norm(x(:));

for k = 1 : size(field, 1)
    f(k,:) = my_dot(ref, field(k,:));
end

