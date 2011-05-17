function [x] = field_update(A, b, C, d)

[x, v, solve_time] = la_quadeq(A'*A, A'*b, C, d);

function [x, v, solve_time] = la_quadeq (A, b, C, d)

n = size(A, 1);
p = size(C, 2);

tic;
w = A \ [C, b];
v = (C' * w(:,1:p)) \ (C' * w(:,end) - d);
x = w(:,end) - w(:,1:p) *  v;
solve_time = toc;
