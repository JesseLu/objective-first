clc; clf; close all;

% Create the initial structure, a silicon waveguide in air.
eps_rows = 40;
eps_cols = 80;
num_eps_pixels = eps_rows * eps_cols;
eps0 = ones([eps_rows eps_cols]);

% Final simulation space size
simulation_width = 160;
simulation_height = 120;

% The maximum and minimum values of epsilon allowed in the structure
max_eps_value = 12.25;
min_eps_value = 1.0;
max_eps_array = 12.25 * ones(size(eps0));
min_eps_array = 1.0 * ones(size(eps0));

% The input mode coupling from and the output mode coupling to
input_mode = 1;
output_mode = 2;
design_frequency = 0.15;

% Run the rectangular waveguide completely through the middle of the
% coupler to begin with
eps0(:,30:50) = max_eps_value;

% Fill everything, except for a 2-cell boundary layer, with epsilon = 9.
% The 2-cell boundary layer is left unchanged because the boundary conditions
% use these cells.
eps0(3:end-2,3:end-2) = 9.0;

% Setup the design. Specifically, we want
% *   a design frequency of 0.15,
% *   to limit epsilon from 1 to 12.25, and
% *   to use the fundamental mode as an input on the left and
%     the second-order mode as the output mode on the right.
spec = setup(design_frequency, eps0, [min_eps_value max_eps_value], [input_mode output_mode]);

% Set the initial values of epsilon in the structure to the initial values.
% These will evolve in time over the course of the optimization.
eps = eps0;

upper_eps = max_eps_value * ones(num_eps_pixels, 1);
lower_eps = min_eps_value * ones(num_eps_pixels, 1);

% The desired output for Hz
c = spec.out.Hz;
% Compute the objective function for the optimization.  We are maximizing
% the overlap of the output mode from the structure with the desired output
% waveguide mode specified.  As a minimization problem, this is equivalent
% ot minimizing the negative of that same overlap integral.
compute_objective = @(input_Hz) (-(c'*input_Hz)*conj(c'*input_Hz));
% Partial derivative of the objective function, g, with respect to the
% field variable input_Hz
g_x_partial = @(input_Hz, S_filter) -c'*S_filter*conj(c'*(S_filter*input_Hz));

% Setup gradient descent parameters.  Run for a fixed number of iterations
% with a fixed step size.
step_size = 20;
num_iter = 40;

for n = 1 : 1 : num_iter
    tic;
    [Ex_adj, Ey_adj, Hz_adj, gradient_adj] = ob1_fdfd_adj(spec.omega, eps, spec.in, g_x_partial, spec.bc);
    adjoint_time = toc;

    fprintf('Iteration #%d, Objective Function Value = %f, Computation Time (sec) = %f!\n', ...
        n, compute_objective(Hz_adj(end,:).'), adjoint_time);
       
    eps(3:end-2,3:end-2) = eps(3:end-2,3:end-2) - step_size * (gradient_adj(3:end-2,3:end-2));
    eps(3:end-2,3:end-2) = min(eps(3:end-2,3:end-2), max_eps_array(3:end-2,3:end-2));
    eps(3:end-2,3:end-2) = max(eps(3:end-2,3:end-2), min_eps_array(3:end-2,3:end-2));
end
    
% Simulate the result and analyze the success of the design process.
simulate(spec, eps, [simulation_width simulation_height]);

save_epsilon_filename = sprintf('adjoint_eps_%d_%d.csv', input_mode, output_mode);
csvwrite(save_epsilon_filename, eps);

