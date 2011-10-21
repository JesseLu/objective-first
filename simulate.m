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
%     P_OUT: Non-negative scalar.
%         The power in the desired ouput mode. The input mode is excited
%         with power ~ 1.0. 
%         
%         For accurate efficiency calculations, measure the output power for
%         an unbroken version of the input waveguide, this is the true
%         (accurate) amount of power excited in the input mode.

% Hard-coded parameters.
t_pml = 10; % Thickness of PML.

    
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


[Ex, Ey, Hz] = ob1_fdfd(spec, eps, t_pml, pad); % Solve.

[P_out] = ob1_calc_power(spec, Ex, Ey, Hz, t_pml, pad);
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
