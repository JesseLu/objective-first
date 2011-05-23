function [mode] = mode_solve(eps, omega, dir)


% Make z the propagation the direction.
switch (dir)
    case {'-x', '+x'}
        eps = struct('x', eps.y(:), 'z', eps.x(:));
    case {'-y', '+y'}
        eps = struct('x', eps.x(:), 'z', eps.z(:));
end

global D

% Shortcut to form a derivative matrix.
S = @(sx, sz) shift_mirror(size(eps.x), -[sx sz]); % Mirror boundary conditions.


    %
    % Build the physics matrix. We assume for now that propagation is in the 
    % z-direction.
    %

div = S(0,0) - S(-1,0); % Divergence.
grad = S(1,0) - S(0,0); % Gradient.

% Matrix that define the eigenproblem.
A = grad * D(1./eps.z) * div * D(eps.x) + omega^2 * D(eps.x);

Ex2Ez = @(Ex, beta) 1 ./ (beta * eps.z) .* (div * (eps.x .* Ex));


    %
    % Solve for the waveguide mode.
    %

% Solve the eigenproblem.
[V, D] = eig(full(A));

% Obtain the wave-vector (beta) and the transverse field.
[beta2, ind] = max(diag(D));
beta = sqrt(beta2);
Ex = V(:,ind);

% Simple scheme to remove arbitrary coefficients/phases in the field.
if (abs(max(Ex)) < abs(min(Ex)))
    Ex = -Ex;
else 
    Ex = Ex;
end

% Obtain the longitudinal field.
Ez = Ex2Ez(Ex, beta);


plot([eps.x(:), eps.z(:)], '.-'); % Plot the cut-away structure.

plot([Ex, Ez], '.-');

