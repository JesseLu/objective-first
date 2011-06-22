function [A, b, B, d] = em_physics1(omega)


    %
    % Helper functions for building matrices.
    %

global S_ DIMS_ D_
N = prod(DIMS_); 

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))];  
Hcurl = [   (S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; 

A_spread = 0.5 * [S_(0,0)+S_(1,0); S_(0,0)+S_(0,1)];
A = @(p) [Ecurl, -i*omega*speye(N); i*omega*D_(A_spread*p), Hcurl];


B = @(x) i * omega * D_(x(1:2*N)) * A_spread; 
d = @(x) -Hcurl * x(2*N+1:3*N);

% Physics residual.
f = @(v) 0.5 * norm(A(v.p)*v.x)^2;

% Gradient.
g = @(v) struct('x', A(v.p)'*(A(v.p)*v.x), ...
                'p', (A_spread'*A_spread) \ (B(v.x)'*(B(v.x)*v.p - d(v.x))));

