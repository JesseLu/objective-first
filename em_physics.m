function [f, g] = em_physics(omega, dims)


    %
    % Helper functions for building matrices.
    %

% Shortcut to form a derivative matrix.
s = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(s(0,1)-s(0,0)),  (s(1,0)-s(0,0))];  
Hcurl = [   (s(0,0)-s(0,-1)); -(s(0,0)-s(-1,0))]; 


% Physics residual.
f = @(v) 0.5 * (norm(Ecurl * v.E - i * omega * v.H)^2 + ...
                norm(Hcurl * v.H + i * omega * (v.eps .* v.E))^2);

% Gradient.
g = @(v) struct( ...
    'E', Ecurl' * (Ecurl * v.E - i * omega * v.H), ...  
    'H', Hcurl' * (Hcurl * v.H + i * omega * (v.eps .* v.E)), ...  
    'eps',  conj(i * omega * v.E) .* ...  
        (Hcurl' * (Hcurl * v.H + i * omega * (v.eps .* v.E))));

return
% Primary physics matrix, electromagnetic wave equation.
A = @(p) [Ecurl, -i*omega*speye(N); i*omega*D_(A_spread * p(:)), Hcurl];

B = @(x) i * omega * D_(x(1:2*N)) * A_spread; 
d = @(x) -Hcurl * x(2*N+1:3*N);
