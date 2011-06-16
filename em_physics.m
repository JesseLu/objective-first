function [f, g] = em_physics(dims, eps)


    %
    % Helper functions for building matrices.
    %

global S_ 

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))];  
Hcurl = [   (S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; 

e = [eps.x(:); eps.y(:)];

% Physics residual.
f = @(v) 0.5 * (norm(Ecurl * v.E - i * omega * v.H)^2 + ...
                norm(Hcurl * v.H + i * omega * (e .* v.E))^2);

% Gradient.
g = @(v) struct( ...
    'E', Ecurl' * (Ecurl * v.E - i * omega * v.H), ...  
    'H', Hcurl' * (Hcurl * v.H + i * omega * (e .* v.E)));

return
% Primary physics matrix, electromagnetic wave equation.
A = @(p) [Ecurl, -i*omega*speye(N); i*omega*D_(A_spread * p(:)), Hcurl];

B = @(x) i * omega * D_(x(1:2*N)) * A_spread; 
d = @(x) -Hcurl * x(2*N+1:3*N);
