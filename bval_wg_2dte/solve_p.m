function [p] = solve_p(B, d, p2, eta)

    %
    % Transform the problem to constrain the solution to the real values.
    %

A = B' * B;
b = real(B' * d);


    %
    % Only allow p to vary where eta < 0.
    %

ind = {find(eta < 0), find(eta >= 0)};
for k = 1 : 2
    n = length(ind{k});
    S{k} = sparse(ind{k}, 1:n, ones(n,1), size(A,2), n);
end

% Formulate the equivalent problem.
Ahat = A * S{1};
bhat = b - A * S{2} * p2(ind{2});


    % 
    % Solve for p1 (the values of p allowed to change). And back-out p.
    %

p1 = Ahat \ bhat;

% Final solution.
p = S{1} * p1 + S{2} * p2(ind{2});
