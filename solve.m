function [eps] = solve(spec, max_iters, min_grad, varargin)
% EPS = SOLVE(SPEC, MAX_ITERS, MIN_GRAD, [METHOD])
%
% Description
%     Perform an objective-first optimization for a nanophotonic waveguide
%     coupler.
% 
% Inputs
%     SPEC: Structure.
%         The output of SETUP(). SPEC defines the waveguide coupler problem.
% 
%     MAX_ITERS: Non-negative integer.
%         Maximum number of iterations for which to run the optimization.
% 
%     MIN_GRAD: Non-negative scalar.
%         Minimum value for the norm of the gradient. The optimization will
%         stop once the norm of the gradient is below MIN_GRAD.
% 
%     METHOD: Character string.
%         The name of the numerical method to use to solve the problem.
%         Currently only the 'sep_convex' option is supported.
% 
% Outputs
%     EPS: 2d array.
%         The permittivity values of the final design.
    

