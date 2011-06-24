function [f, g] = em_physics(omega, x)


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
f = @(p) 0.5 * norm(field_template .* (A(p)*x))^2;

% Gradient.
g = @(p) B(x)'*(B(x)*p - d(x));

