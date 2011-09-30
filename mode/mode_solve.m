function [mode] = mode_solve(eps, omega, dir, mode_num)


    %
    % Package/re-align so that z is the propagation direction.
    %

% Make z the propagation the direction.
switch (dir)
    case {'x-', 'x+'}
        eps = struct('x', eps.y(:), 'z', eps.x(:));
    case {'y-', 'y+'}
        eps = struct('x', eps.x(:), 'z', eps.y(:));
end


    %
    % Build the physics matrix. We assume for now that propagation is in the 
    % z-direction.
    %

D_ = @(x) spdiags(x(:), 0, numel(x), numel(x)); % Diagonalization function.

% Shortcut to form a derivative matrix.
S = @(sx, sz) shift_circle(size(eps.x), -[sx sz]); % Mirror boundary conditions.

div = S(0,0) - S(-1,0); % Divergence.
grad = S(1,0) - S(0,0); % Gradient.

% Matrix that defines the eigenproblem.
A = grad * D_(1./eps.z) * div * D_(eps.x) + omega^2 * D_(eps.x);

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
[beta2, ind] = sort(real(diag(C)), 'descend');
beta = sqrt(beta2(mode_num));
Ex = V(:,ind(mode_num));

% plot([eps.x(:), eps.z(:)], '.-'); pause
% for k = 1 : size(V, 2)
%     plot(V(:,k))
%     fprintf('%d: %f\n', k, (beta2(k)))
%     pause
% end

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


    %
    % Normalize power (Poynting vector) to 1.
    %

P = abs(sum(sign(eps.x) .* Ex .* Hy));
Ex = Ex ./ sqrt(P);
Ez = Ez ./ sqrt(P);
Hy = Hy ./ sqrt(P);


    %
    % Re-package for the actual propagation direction.
    %

mode = struct('beta', beta, 'Et', Ex, 'El', Ez, 'Ht', Hy);
