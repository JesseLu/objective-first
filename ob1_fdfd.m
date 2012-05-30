function [Ex, Ey, Hz] = ob1_fdfd(omega, eps, in, bc)
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

t_pml = 20;

% Expand eps to include room for the pml padding
if strcmp(bc, 'pml') 
    eps = ob1_pad_eps(eps, t_pml * [1 1 1 1]);
elseif strcmp(bc, 'per') 
    eps = ob1_pad_eps(eps, t_pml * [1 1 0 0]);
end
dims = size(eps); % New size of the structure.
[eps_x, eps_y] = ob1_interp_eps(eps); % Obtain x- and y- components of eps.

    
    %
    % Form the system matrix which will be solved.
    %
    
    %
    % Determine the input excitation.
    %

b = zeros(dims); % Input excitation, equivalent to magnetic current source.
in_pos = t_pml+2; % Cannot be 1, because eps interpolation wreaks havoc at border.

% For one-way excitation in the forward (to the right) direction,
% we simple cancel out excitation in the backward (left) direction.
if strcmp(bc, 'pml')
    b_pad = [floor((dims(2) - length(in.Hz))/2), ...
            ceil((dims(2) - length(in.Hz))/2)];
elseif strcmp(bc, 'per')
    b_pad = [0 0];
end
b(in_pos+1, b_pad(1)+1:end-b_pad(2)) = in.Hz;
b(in_pos, b_pad(1)+1:end-b_pad(2)) = -in.Hz * exp(i * in.beta);

b = b ./ eps_y; % Convert from field to current source.

% Solve.
[Ex, Ey, Hz] = ob1_fdfd_adv(omega, {eps_x, eps_y}, ones(dims), reshape(b, dims), bc, t_pml);

% Cut off the pml parts.
if strcmp(bc, 'pml')
    Ex = Ex(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
    Ey = Ey(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
    Hz = Hz(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
elseif strcmp(bc, 'per')
    Ex = Ex(t_pml+1:end-t_pml, :);
    Ey = Ey(t_pml+1:end-t_pml, :);
    Hz = Hz(t_pml+1:end-t_pml, :);
end
