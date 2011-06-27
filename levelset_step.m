function [phi] = levelset_step(phi, dphi, s)
% PHI = LEVELSET_STEP(PHI, DPHI, S)
% 
% Description
%     Update the interface with the suggested step size s. Take into account 
%     the numerical stability (CFL) condition as well.

% Compute the Hamiltonian for motion in the normal direction.
H_normal = dphi .* norm_gradient(phi, dphi);

nb = @(x) x(:) .* ((phi(:) >= -3) & (phi(:) <= 3));
s_max = 1 / max(abs(nb(H_normal)));

% Update phi.
phi = phi - min([s, s_max]) * H_normal;
% fprintf('phi step: %e\n', min([s, s_max]));

% [phi, err] = signed_distance(phi, 1e-2); % Make phi more sdf-like.
