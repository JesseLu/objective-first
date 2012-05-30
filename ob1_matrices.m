function [A, S] = ob1_matrices(dims)
% Form matrices used for objective-first method.

    %
    % Form A matrices.
    %

% Shortcut to form a 2D derivative matrix.
S_ = @(sx, sy) ob1_shift_matrix(dims, -[sx sy]); % Mirror boundary conditions.

A{1} = [-(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))]; % Curl for E-field.
A{2} = 0.5 * [S_(0,0)+S_(1,0); S_(0,0)+S_(0,1)]; % Spread for epsilon.
A{3} = [(S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; % Curl for H-field.

    
    %
    % Form selection matrices, used to select either interior,
    % or boundary elements of x and p.
    %

[S.int, S.bnd] = my_selection(dims, [2 2]); % 2-cell border on all four sides.
[S.res] = my_selection(dims, [1 1]); % Used to calculate the residual.



function [int, bnd] = my_selection(dims, border)
% Form selection matrices for both interior and boundary values.

% Temporary field which we use to select boundary elements.
field = zeros(dims);
field(border(1)+1:end-border(1), border(2)+1:end-border(2)) = 1;

% Find which elements are either on the interior or on the boundary.
ind_int = find(field == 1);
ind_bnd = find(field == 0);

% Function which forms the sparse matrix.
my_sparse = @(ind) sparse(1:length(ind), ind, ones(length(ind),1), ...
                            length(ind), prod(dims));

% For the selection matrices.
int = my_sparse(ind_int);
bnd = my_sparse(ind_bnd);
