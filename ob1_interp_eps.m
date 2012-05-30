function [eps_x, eps_y] = ob1_interp_eps(eps)
% Interpolate eps.
%
% We actually interpolate the inverse of eps since this is what occurs in
% the SOLVE() function.
%
% Lastly, we wrap-around (assume periodic boundary conditions).

eps_x = (0.5 * (eps.^-1 + eps([2:end, 1], :).^-1)).^-1;
eps_y = (0.5 * (eps.^-1 + eps(:, [2:end, 1]).^-1)).^-1;
