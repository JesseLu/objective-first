function [Ahat, bhat, add_border] = setup_border_insert(A, x)

% Helper variables.
global DIMS_
dims = DIMS_;
N = prod(dims);


    %
    % Find the border and interior indices.
    %

% Create a temporary array with ones on the border and zeros on the interior.
temp_field = zeros(dims);
temp_field(:,[1 end]) = 1;
temp_field([1 end],:) = 1;
temp_field = repmat(temp_field, 1, 3);

% Find the border and the interior indices.
bor = find(temp_field(:) == 1);
int = find(temp_field(:) == 0);


    %
    % Truncate the physics matrix, and create the vector b.
    %

Ahat = A(:,int);
bhat = A(:,bor) * x(bor);


    %
    % Simple function to re-insert the border values.
    %

n_bor = length(bor);
n_int = length(int);
add_border = @(xhat) ...
    sparse(bor, 1:n_bor, ones(n_bor, 1), 3*N, n_bor) * x(bor) + ...
    sparse(int, 1:n_int, ones(n_int, 1), 3*N, n_int) * xhat;
    

