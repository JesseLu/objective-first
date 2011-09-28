function [mode] = mode_solve(eps, omega, dir)


    %
    % Package/re-align so that z is the propagation direction.
    %

% Make z the propagation the direction.
switch (dir)
    case {'-x', '+x'}
        eps = struct('x', eps.y(:), 'z', eps.x(:));
    case {'-y', '+y'}
        eps = struct('x', eps.x(:), 'z', eps.y(:));
end


    %
    % Build the physics matrix. We assume for now that propagation is in the 
    % z-direction.
    %

global D % Diagonalization function.

% Shortcut to form a derivative matrix.
S = @(sx, sz) shift_mirror(size(eps.x), -[sx sz]); % Mirror boundary conditions.

div = S(0,0) - S(-1,0); % Divergence.
grad = S(1,0) - S(0,0); % Gradient.

% Matrix that defines the eigenproblem.
A = grad * D(1./eps.z) * div * D(eps.x) + omega^2 * D(eps.x);

% Function to obtain the longitudinal component.
Ex2Ez = @(Ex, beta) i ./ (beta * eps.z) .* (div * (eps.x .* Ex));

% Function to obtain the Hy-field.
E2Hy = @(Ex, Ey, beta) 1 / (i*omega) * (i * beta * Ex - grad * Ey);


    %
    % Solve for the waveguide mode.
    %

% Solve the eigenproblem.
[V, C] = eig(full(A));

% Obtain the wave-vector (beta) and the transverse field.
[beta2, ind] = max(diag(C));
beta = sqrt(beta2);
Ex = V(:,ind);

% Simple scheme to remove arbitrary coefficients/phases in the field.
if (abs(max(Ex)) < abs(min(Ex)))
    Ex = -Ex;
else 
    Ex = Ex;
end

% Obtain the longitudinal E-field (Ez).
Ez = Ex2Ez(Ex, beta);

% Obtain the transverse H-field (Hy).
Hy = E2Hy(Ex, Ez, beta);

% Info for sourcing.
Jx = Ex .* eps.x; 
Jz = Ez .* eps.z; 
My = Hy;


    %
    % Re-package for the actual propagation direction.
    %

mode = struct('beta', beta,     'Et', Ex, 'El', Ez, 'Ht', Hy, ...
                                'Jt', Jx, 'Jl', Jz, 'Mt', My);
