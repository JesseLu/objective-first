function [x] = arrow_solve(d, A, U, V, b)
% Solve a linear system which consists of the sum of two matrices:
%     ([diag(d) A^T; A 0] + U*V^T) * x = b.


% Rank of second matrix (assume columns are linearly independent).
m = size(U, 2); 


    %
    % Calculate elements needed to use the matrix inversion lemma.
    %

% Calculate solution to the first (arrow matrix)
w = my_block_elim(d, A, [U, b]); 

% Separate out the two components needed to use the lemma.
Y = w(:,1:m); % Solution of A_hat^-1 U.
z = w(:,m+1:end); % Solution of A_hat^-1 b.

    
    % 
    % Apply matrix inversion lemma to obtain solution.
    %

x = z - Y * inv(eye(m) + V' * Y) * (V' * z);


function [x] = my_block_elim(d, A, b)
% Uses block elimination to solve an arrow matrix.
% Assumes that A is skinny and with linearly independent columns.

% Separate out upper and lower blocks of b.
b1 = b(1:size(A,2),:);
b2 = b(size(A,2)+1:end,:);

D_inv = spdiags(d.^-1, 0, length(d), length(d)); % Inverse of diagonal block.
S_inv = -inv(full(A * D_inv * A')); % Inverse Schur complement (1,1 block).

x2 = S_inv * (b2 - A * D_inv * b1); % Solve for lower block.
x1 = D_inv * (b1 - A' * x2); % Solve for upper block.

x = [x1; x2]; % The solution.

