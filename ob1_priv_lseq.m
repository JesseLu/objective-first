function [x] = ob1_priv_lseq(A, b, C, d, eta)

Ahat = [A'*A+eta*speye(size(A,2)), C'; C, zeros(size(C,1))];
bhat = [A'*b; d];

z = Ahat \ bhat;
x = z(1:size(A,2));
v = z(size(A,2)+1:end);
