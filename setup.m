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
%         The initial permittivity of the structure.
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
%     PHASE: Scalar (optional).
%         The relative phase of the output mode relative to the input mode.
%         Default value assumes that the input and output modes each propagate 
%         through half of the design space.
% 
% Outputs
%     SPEC: Structure.
%         Contains the input and output information which defines the design 
%         problem. The input and output waveguide modes have been normalized 
%         to have a Poynting vector (power throughput) of unity.
%         SPEC is mainly used as an input to the DESIGN() function.
% 
% Examples

dims = size(eps0);

% Create the specification structure.
spec.omega = omega;
spec.eps0 = eps0;
spec.in.eps = eps0(1,:);
spec.out.eps = eps0(end,:);
spec.eps_lims = sort(eps_lims);


    %
    % Calculate the input and output modes.
    %

figure(1);
fprintf('Input mode calculation (figure 1)\n');
[spec.in.beta, spec.in.Hz, spec.in.Ey] = ...
    ob1_wgmode(spec.omega, spec.in.eps, mode_nums(1));

figure(2);
fprintf('Output mode calculation (figure 2)\n');
[spec.out.beta, spec.out.Hz, spec.out.Ey] = ...
    ob1_wgmode(spec.omega, spec.out.eps, mode_nums(2));



