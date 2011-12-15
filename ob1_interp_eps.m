function [eps_x, eps_y, A_spread] = ob1_interp_eps(eps)
% Interpolate eps.
%
% We actually interpolate the inverse of eps since this is what occurs in
% the SOLVE() function.
%
% Lastly, we wrap-around (assume periodic boundary conditions).

eps_x = (0.5 * (eps + eps([2:end, 1], :)));
eps_y = (0.5 * (eps + eps(:, [2:end, 1])));

% Shortcut to form a 2D derivative matrix.
S_ = @(sx, sy) ob1_shift_matrix(size(eps), -[sx sy]); % Mirror boundary conditions.

A_spread = 0.5 * [(S_(0,0)+S_(1,0)); (S_(0,0)+S_(0,1))];
