function [mode] = find_wg_mode(p, omega, device, side)

    %
    % Find the relevant slice of the structure.
    %

% Device boundaries.
lim.x = device.offset(1) + [0 device.dims(1)-1];
lim.y = device.offset(2) + [0 device.dims(2)-1];

% "Cut out" the structure.
switch (side)
    case '-x'
        p = p(  lim.x(1),           lim.y(1):lim.y(2));
    case '+x'
        p = p(  lim.x(2),           lim.y(1):lim.y(2));
    case '-y'
        p = p(  lim.x(1):lim.x(2),  lim.y(1));
    case '+y'
        p = p(  lim.x(1):lim.x(2),  lim.y(2));
end


    %
    % Build the physics matrix. We assume for now that propagation is in the 
    % z-direction.
    %

% Helper function to cleanly define derivatives.
s = @(sx, sy) mirror_shift(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
my_diag = @(x) spdiags(x(:), 0, numel(x), numel(x));

div = s(0,0) - s(-1,0); % Divergence.
grad = s(1,0) - s(0,0); % Gradient.

A = grad * my_diag(

