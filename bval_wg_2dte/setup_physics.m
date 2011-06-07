function [A, B, d] = setup_physics(omega, p2e)


    %
    % Helper functions for building matrices.
    %

% Allows for a clean definition of curl.
global S_ D_ DIMS_

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))];  
Hcurl = [   (S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; 

% Find the indices that will compose the border.
N = prod(DIMS_);



    % 
    % Build the matrices that we will be using.
    %

% Primary physics matrix, electromagnetic wave equation.
A = @(p) [Ecurl, -i*omega*speye(N); i*omega*D_(p2e(p)), Hcurl];

B = @(x) i * omega * D_(x(1:2*N)); 
B = @(x) i * omega * (x(1:2*N)); 
d = @(x) -Hcurl * x(2*N+1:3*N);
