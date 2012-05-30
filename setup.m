function [spec] = setup(omega, eps0, eps_lims, mode_nums, varargin)
% SPEC = SETUP(OMEGA, EPS0, EPS_LIMS, MODE_NUMS, [PHASE])
%
% Description
%     Sets up the nanophotonic waveguide coupler design problem.
% 
%     Specifically, SETUP() determines the boundary values at the 
%     input and output edges of the device (left and right respectively).
% 
% Inputs
%     OMEGA: Positive scalar.
%         The design (angular) frequency. 
% 
%     EPS0: 2-d array.
%         The initial permittivity of the structure. It is recommended that
%         the permittivities in the two leftmost and two rightmost layers of
%         cells be uniform, since the fields are fixed in those cells.
%     
%     EPS_LIMS: 2-element vector.
%         The minimum and maximum allowable values for the permittivity.
% 
%     MODE_NUMS: 2-element vector of integers.
%         The order of the input and output waveguide modes, respectively.
%         For example, to select the fundamental mode for the input and
%         the second-order mode for the output use MODE_NUMS = [1 2].
% 
%         To interactively select a mode, input a mode number of 0.
% 
%         If the selected modes are non-propagating (evanesent), a warning will 
%         be issued.
% 
%     BOUNDARY: String (optional).
%         If set to 'periodic', then we assume the top and bottom boundaries to 
%         to have periodic boundary conditions.
% 
% Outputs
%     SPEC: Structure.
%         Contains the input and output information which defines the design 
%         problem. Sets up the boundary-value formulation for use as the 
%         design objective. 
% 
%         SPEC is mainly used as an input to the DESIGN() function.

dims = size(eps0);

% Create the specification structure.
spec.omega = omega;
spec.eps0 = eps0;
spec.eps_lims = sort(eps_lims);


    %
    % Calculate the input and output modes.
    %

figure(1);
fprintf('Input mode calculation (figure 1)\n');
[spec.in.beta, spec.in.Hz, spec.in.Ey] = ...
    ob1_wgmode(spec.omega, eps0(1,:), mode_nums(1));

figure(2);
fprintf('Output mode calculation (figure 2)\n');
[spec.out.beta, spec.out.Hz, spec.out.Ey] = ...
    ob1_wgmode(spec.omega, eps0(end,:), mode_nums(2));


    %
    % Set up the boundary values.
    %

% Determine phase relation between input and output modes.
phase = mean([spec.in.beta, spec.out.beta]) * (size(eps, 1) - 1);
if isempty(varargin)
    spec.bc = 'pml';
elseif strcmp(varargin{1}, 'periodic')
    spec.bc = 'per';
else
    error('Not a valid option for BOUNDARY.');
end

% Create boundary field conditions.
% Two layers of Hz fields are needed since we must fix Hz as well as its 
% derivative in the normal direction.
spec.Hz0 = zeros(dims);
spec.Hz0(1,:) = spec.in.Hz;
spec.Hz0(2,:) = spec.in.Hz * exp(i * spec.in.beta); 
spec.Hz0(end-1,:) = spec.out.Hz * exp(i * (phase - spec.out.beta)); 
spec.Hz0(end,:) = spec.out.Hz * exp(i * phase); 
