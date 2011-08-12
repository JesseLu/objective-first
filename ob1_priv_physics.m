function [A, B, d, A_spread] = ob1_priv_physics(omega, dims)
% 
% Description
%     Form the underlying matrices (and vectors) that represent the physics.

N = prod(dims);

% Shortcut to form a derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
D_ = @(x) spdiags(x(:), 0, numel(x), numel(x));

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))];  
Hcurl = [   (S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; 

% Matrix that determines x- and y-components of epsilon from the z-component.
A_spread = 0.5 * [S_(0,0)+S_(1,0); S_(0,0)+S_(0,1)];

% Matrix for the field.
% A = @(epsilon) [Ecurl, -i*omega*speye(N); i*omega*D_(A_spread*epsilon), Hcurl];
A = @(epsilon) Ecurl * D_(A_spread*epsilon) * Hcurl - omega^2 * speye(N);

% Matrix and vector for structure.
B = @(x) Ecurl * D_(Hcurl * x) * A_spread; 
d = @(x) omega^2 * x;

