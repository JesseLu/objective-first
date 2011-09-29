function [H] = sim(omega0, epsilon)

dims = size(epsilon);
N = prod(dims);

% Shortcut to form a 2D derivative matrix.
S_ = @(sx, sy) shift_circle(dims, -[sx sy]); % Mirror boundary conditions.

A{1} = [-(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))]; % Curl for E-field.
A{2} = 0.5 * [S_(0,0)+S_(1,0); S_(0,0)+S_(0,1)]; % Spread for epsilon.
A{3} = [(S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; % Curl for H-field.

% Helper function for sparse diagonal matrix.
D_ = @(z) sparse(1:length(z), 1:length(z), z, length(z), length(z));

p = epsilon(:).^-1;

A = A{1} * D_(A{2} * p) * A{3};

[V, D] = eig(full(A));

d = sqrt(diag(D));
[temp, ind] = sort(abs(real(d) - omega0));

for k = 1 : length(d)
    imagesc(reshape(V(:,ind(k)), dims)'); axis equal tight
    d(ind(k))
    pause
end


function my_plot(z)

c = max(abs(z(:)));
imagesc(z', c * [-1 1]); axis equal tight
