function [eps, Ex, Ey, Hz] = ob1_fdfd(omega, eps, in)
% 
% Description
%     Solve a FDFD (finite-difference, frequency-domain) problem, using the
%     input mode as the source term.
% 
%     OB1_FDFD actually expands the structure in order to add absorbing
%     boundary layers (pml), simulates the structure by sourcing it from
%     within the pml, and then truncates the solution so the user does not
%     see the absorbing pml region.
% 
%     Note that a time dependence of exp(-i * omega * t) is assumed for all
%     fields.



[A, b, Ecurl, eps] = get_Ab(omega, eps, in);

    %
    % Solve.
    %

% % Test to make sure get_Ab and get_Bd work.
% E = rand(2*numel(eps), 1);
% [B, d] = get_Bd(omega, E, size(eps), b);
% rx = A * E - b;
% rp = B * eps(:) - d;
% norm(rx - rp)
% norm(rx)
% norm(rp)

% Hz = A \ b; % This should be using sparse matrix factorization. 
E = A \ b; % This should be using sparse matrix factorization. 


Hz = -(i/omega) * Ecurl * E; % Obtain Hz-field.

% Reshape and extract all three fields.
dims = size(eps);
Ex = reshape(E(1:prod(dims)), dims);
Ey = reshape(E(prod(dims)+1:end), dims);
Hz = reshape(Hz, dims);

% % Cut off the pml parts.
% Ex = Ex(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
% Ey = Ey(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
% Hz = Hz(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
