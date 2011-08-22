function [x] = ob1_priv_lseq(A, b, C, d)

Ahat = [A'*A, C'; C, zeros(size(C,1))];
bhat = [A'*b; d];

z = Ahat \ bhat;
x = z(1:size(A,2));
v = z(size(A,2)+1:end);
