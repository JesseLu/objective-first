function [A, S] = ob1_setup_matrices(dims, active_box)

path(path, '~/wave-tools/helper');

    %
    % Form A matrices.
    %

% Shortcut to form a 2D derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

A{1} = [-(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))]; % Curl for E-field.
A{2} = 0.5 * [S_(0,0)+S_(1,0); S_(0,0)+S_(0,1)]; % Spread for epsilon.
A{3} = [(S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; % Curl for H-field.

    
    %
    % Form selection matrices, used to select active elements of x and p.
    %

my_diag = @(x) sparse(1:numel(x), 1:numel(x), x(:), numel(x), numel(x));

S.x = zeros(dims);
S.x(2:end-1, 2:end-1) = 1;
S.x = my_diag(S.x);

b = ceil((dims - active_box) / 2);
S.p = zeros(dims);
S.p(b(1)+1:end-b(1), b(2)+1:end-b(2)) = 1;
S.p = my_diag(S.p);
