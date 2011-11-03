function [x] = arrow_solve(d, A, U, V, b)
% Solve a linear system which consists of the sum of two matrices:
%     ([diag(d) A^T; A 0] + U*V^T) * x = b.

m = size(U, 2);

w = my_block_elim(d, A, [U, b]);
Y = w(:,1:m); % Solution of A_hat^-1 U.
z = w(:,m+1:end); % Solution of A_hat^-1 b.

x = z - Y * inv(eye(m) + V' * Y) * (V' * z); % Use matrix inversion lemma.

function [x] = my_block_elim(d, A, b)
b1 = b(1:size(A,2),:);
b2 = b(size(A,2)+1:end,:);

D_inv = spdiags(d.^-1, 0, length(d), length(d));

S_inv = -inv(full(A * D_inv * A')); % Inverse Schur complement (1,1 block).

% Solve for lower block.
x2 = S_inv * (b2 - A * D_inv * b1); 
x1 = D_inv * (b1 - A' * x2);

x = [x1; x2];

