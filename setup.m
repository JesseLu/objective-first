function [spec] = setup(omega, eps_in, order_in, eps_out, order_out)
% SPEC = SETUP(OMEGA, EPS_IN, ORDER_IN, EPS_OUT, ORDER_OUT)
%
% Description
%     Sets up the nanophotonic waveguide coupler design problem.
% 
%     Specifically, SETUP() determines the input and output waveguide modes,
%     based on the dielectric structure of the input and output waveguides.
% 
% Inputs
%     OMEGA: Positive scalar.
%         The design (angular) frequency. 
% 
%     EPS_IN, EPS_OUT: 1-d arrays.
%         The (relative) permittivities for the input and output waveguides,
%         respectively.
% 
%     ORDER_IN, ORDER_OUT: Positive integers.
%         The order of the input and output waveguide modes, respectively.
%         For example, to select the fundamental mode, use ORDER_IN or 
%         ORDER_OUT = 1. For second-order mode, use ORDER_IN or ORDER_OUT = 2.
%         To interactively select modes, user ORDER_IN or ORDER_OUT = [].
% 
%         If the selected mode is non-propagating (evanesent), a warning will 
%         be issued.
% 
% Outputs
%     SPEC: Structure.
%         Contains the input and output information which defines the design 
%         problem. The input and output waveguide modes have been normalized 
%         to have a Poynting vector (power throughput) of unity.
%         SPEC is mainly used as an input to the DESIGN() function.
% 
% Examples

% Create the specification structure.
spec.omega = omega;
spec.in.eps = reshape(eps_in, 1, length(eps_in));
spec.out.eps = reshape(eps_out, 1, length(eps_out));


fprintf('Input mode calculation\n');
[spec.in.beta, spec.in.Hz, spec.in.Ey] = ...
    ob1_wgmode(spec.omega, spec.in.eps, order_in);

fprintf('Output mode calculation\n');
[spec.out.beta, spec.out.Hz, spec.out.Ey] = ...
    ob1_wgmode(spec.omega, spec.out.eps, order_out);



