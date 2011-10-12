function [P_out] = simulate(spec, eps, dims)
% P_OUT = SIMULATE(SPEC, EPS, DIMS)
% 
% Description
%     Simulates the design and determines the performance.
% 
% Inputs
%     SPEC: Structure.
%         This is the specification of the problem, as obtained from SETUP().
% 
%     EPS: 2-d array.
%         The permittivity values of the design to be simulated.
% 
%     DIMS: 2-element vector.
%         The size of the simulation in the x- and y- directions respectively.
%         These values should be considerable larger than size(EPS).
% 
% Outputs
%     EFF: Non-negative scalar.
%         The coupling efficiency where 1.0 is 100% conversion efficiency.
% 
% Examples

% Hard-coded parameters.
t_pml = 10; % Thickness of PML.
sigma_pml = 1 / spec.omega; % Strength of PML.
exp_pml = 2.5; % Exponential spatial increase in pml strength.

    
    % 
    % Determine epsilon for the simulation.
    %

% Number of extra cells to add to eps.
pad = [ floor((dims(1) - size(eps, 1))/2), ...
        ceil((dims(1) - size(eps, 1))/2), ...
        floor((dims(2) - size(eps, 2))/2), ...
        ceil((dims(2) - size(eps, 2))/2)];

% Expand eps to the full simulation size.
eps = cat(1, repmat(eps(1,:), pad(1), 1), eps, repmat(eps(end,:), pad(2), 1));
eps = cat(2, repmat(eps(:,1), 1, pad(3)), eps, repmat(eps(:,end), 1, pad(4)));

[eps_x, eps_y] = ob1_interp_eps(eps); % Get x and y components of eps.


    %
    % Build the simulation matrix.
    %

% Shortcut to form a derivative matrix.
S = @(sx, sy) ob1_shift_matrix(dims, -[sx sy]);

% Helper function to create stretched-coordinate PML absorbing layers.
scx = @(sx, sy) ob1_stretched_coords(dims, [1 dims(1)+0.5], [sx, sy], ...
    'x', t_pml, sigma_pml, exp_pml);
scy = @(sx, sy) ob1_stretched_coords(dims, [1 dims(2)+0.5], [sx, sy], ...
    'y', t_pml, sigma_pml, exp_pml);

% Define the curl operators as applied to E and H, respectively.
Ecurl = [scy(.5,.5)*-(S(0,1)-S(0,0)), scx(.5,.5)*(S(1,0)-S(0,0))];  
Hcurl = [scy(.5,0)*(S(0,0)-S(0,-1));  scx(0,.5)*-(S(0,0)-S(-1,0))]; 

% Diagonal matrix for 1/epsilon.
inv_eps = spdiags([eps_x(:).^-1; eps_y(:).^-1], 0, 2*prod(dims), 2*prod(dims));

% This is the matrix that we will solve.
A = Ecurl * inv_eps * Hcurl - spec.omega^2 * speye(prod(dims));

    
    %
    % Determine the input excitation.
    %

b = zeros(dims); % Input excitation, equivalent to magnetic current source.
in_pos = max([t_pml+1, round(pad(1)/2)]); % Location of input excitation.

% For one-way excitation in the forward (to the right) direction,
% we simple cancel out excitation in the backward (left) direction.
b(in_pos+1, pad(3)+1:end-pad(4)) = spec.in.Hz;
b(in_pos, pad(3)+1:end-pad(4)) = -spec.in.Hz * exp(i * spec.in.beta);

b = b ./ eps_y; % Convert from field to current source.

% Normalization factor so that the input power is unity.
b = -i * 2 * spec.in.beta / (1 - exp(i * 2 * spec.in.beta)) *  b;

b = b(:); % Vectorize.


    %
    % Solve.
    %

Hz = A \ b; % This should be using sparse matrix factorization. 

E = 1/spec.omega * inv_eps * Hcurl * Hz; % Obtain E-fields.

% Reshape and extract all three fields.
Ex = reshape(E(1:prod(dims)), dims);
Ey = reshape(E(prod(dims)+1:end), dims);
Hz = reshape(Hz, dims);


    %
    % Calculate power output to desired mode.
    %
        
% Location for power calculation.
out_pos = min([round(dims(1)-pad(2)/2), dims(1)-t_pml-1]); 

% Project y onto x.
proj = @(x, y) (dot(x(:), y(:)) / norm(x(:))^2) * x(:);

% Calculate the power in the desired output mode.
P_out = abs(dot(proj(spec.out.Hz, Hz(out_pos,pad(3)+1:end-pad(4))), ...
                proj(spec.out.Ey * exp(i * 0.5 * spec.out.beta), ...
                    Ey(out_pos,pad(3)+1:end-pad(4)))));
            

    %
    % Print and plot results.
    %

fprintf('Output power in desired mode (input power approx. 1.0) : %1.3f\n', ...
    P_out);

ob1_plot(dims, {'\epsilon', eps}, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});

% % The following commands may be used (uncommented) in order to plot more
% % field information.
% % figure(1); 
% ob1_plot(dims, {'\epsilon', eps}, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});
% 
% % Plot all fields.
% figure(2); 
% ob1_plot(dims, ...
%     {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, {'Re(Hz)', real(Hz)}, ...
%     {'Im(Ex)', imag(Ex)}, {'Im(Ey)', imag(Ey)}, {'Im(Hz)', imag(Hz)}, ...
%     {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});
% 
% % Plot absolute value of all three fields.
% figure(3);
% ob1_plot(dims, {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});
