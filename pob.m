function [f, g] = pob(A, b, x)
f = 0.5 * norm(A * x - b)^2;
g = A' * (A * x - b);
