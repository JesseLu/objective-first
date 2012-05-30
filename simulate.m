function [eff, eps, Ex, Ey, Hz] = simulate(spec, eps, dims)
% [EFF, EPS, EX, EY, HZ] = SIMULATE(SPEC, EPS, DIMS)
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
%     EFF : Non-negative scalar.
%         Conversion efficiency from input to output mode. EFF is calculated by
%         first exciting an unbroken waveguide and measuring its power output.
%         This number is then taken as the input power for the actual 
%         simulation, where the coupler device is used.
%         
%     EPS, EX, EY, HZ: 2-d arrays.
%         Simulation values of the 2D TE fields. The fields will be of size 
%         DIMS. Perfectly-matched layers (PML) absorbing boundaries are used
%         in the simulation, but the field values are excluded from the 
%         EPS, EX, EY, and HZ output values.

    
    % 
    % Determine epsilon for the simulation.
    %

% Number of extra cells to add to eps.
if strcmp(spec.bc, 'per') % Override for periodic case.
    dims(2) = size(eps, 2);
end
pad = [ floor((dims(1) - size(eps, 1))/2), ...
        ceil((dims(1) - size(eps, 1))/2), ...
        floor((dims(2) - size(eps, 2))/2), ...
        ceil((dims(2) - size(eps, 2))/2)];
x_out = round(dims(1)-pad(2)/2+1) : dims(1);
y_out = pad(3)+1 : dims(2)-pad(4);

eps = ob1_pad_eps(eps, pad);


    %
    % Determine input power by solving a structure with only an input
    % waveguide.
    %

[Ex, Ey, Hz] = ob1_fdfd(spec.omega, repmat(eps(1,:), size(eps,1), 1), spec.in, spec.bc);

ob1_plot(size(eps), {'\epsilon', eps}, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});
P_in = ob1_calc_power(Ey(x_out,y_out), Hz(x_out,y_out), spec.in);


    %
    % Calculate output power for the structure (and actual output mode).
    %

[Ex, Ey, Hz] = ob1_fdfd(spec.omega, eps, spec.in, spec.bc);
[P_out, err, filtered_percent] = ...
    ob1_calc_power(Ey(x_out,y_out), Hz(x_out,y_out), spec.out);


    %
    % Print and plot results.
    %

eff = P_out / P_in;
fprintf('P_in: %1.3e \nP_out: %1.3e \neff: %1.3f%%\n', ...
    P_in, P_out, 100 * eff);

ob1_plot(dims, {'\epsilon', eps}, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});
% ob1_plot(size(Hz), {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});

% The following commands may be used (uncommented) in order to plot more
% field information.

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
