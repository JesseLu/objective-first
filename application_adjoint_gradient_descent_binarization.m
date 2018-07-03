clc; clf; close all;

% Read in an already optimized structure for continuous values of epsilon
% and we will start from there in the binarization.
eps = csvread('adjoint_eps_1_2.csv');
[eps_rows, eps_cols] = size(eps);
num_eps_pixels = eps_rows * eps_cols;

% Final simulation space size
simulation_width = 160;
simulation_height = 120;

max_eps_value = 12.25;
min_eps_value = 1.0;
eps_range = max_eps_value - min_eps_value;
eps_midpt = 0.5 * (max_eps_value + min_eps_value);
max_eps_array = max_eps_value * ones(size(eps));
min_eps_array = min_eps_value * ones(size(eps));

input_mode = 1;
output_mode = 2;
design_frequency = 0.15;

% Setup the design. Specifically, we want
% *   a design frequency of 0.15,
% *   to limit epsilon from 1 to 12.25, and
% *   to use the fundamental mode as an input on the left and
%     the second-order mode as the output mode on the right.
spec = setup(design_frequency, eps0, [min_eps_value max_eps_value], [input_mode output_mode]);

% The desired output for Hz
c = spec.out.Hz;

compute_objective = @(input_Hz) (-(c'*input_Hz)*conj(c'*input_Hz));
g_x_partial = @(input_Hz, S_filter) -c'*S_filter*conj(c'*(S_filter*input_Hz));

rescale = @(x) (x - min_eps_value) / eps_range;
scale = @(x) (min_eps_value + x * eps_range);

top = @(x, beta, eta) (tanh(beta * eta) + tanh(beta * (x - eta)));
bottom = @(x, beta, eta) ones(size(x))*(tanh(beta * eta) + tanh(beta * (1 - eta)));
new_func = @(x, beta, eta) top(x, beta, eta) ./ bottom(x, beta, eta);
apply_new_func = @(x, beta, eta) scale(new_func(rescale(x), beta, rescale(eta)));

new_func_deriv = @(x, beta, eta) beta*power(sech(beta*(x-eta)), 2) ./ bottom(x, beta, eta);
apply_new_func_deriv = @(x, beta, eta) new_func_deriv(rescale(x), beta, eta);

num_iter = 20;
strength = 0.01;
min_strength = 0.01;
max_strength = 20.0;
% Step the strength linearly each iteration
strength_values = linspace(min_strength, max_strength, num_iter);

step_size = 20;

for n = 1 : 1 : num_iter
    strength = strength_values(n);

    tic;
    z = reshape(apply_new_func(eps(:), strength, eps_midpt), size(eps));

    [Ex_adj, Ey_adj, Hz_adj, gradient_adj] = ob1_fdfd_adj(spec.omega, z, spec.in, g_x_partial, spec.bc);
    
    dz_de = zeros(eps_rows * eps_cols, 1);
    flatten_eps = eps(:);
    for m = 1 : 1 : num_eps_pixels
        dz_de(m) = apply_new_func_deriv(flatten_eps(m), strength, center);
    end
  
    gradient_adj = reshape(diag(dz_de) * gradient_adj(:), size(eps));
        
    adjoint_time = toc;

    cur_objective = compute_objective(Hz_adj(end,:).');

    fprintf('Iteration #%d, Sigmoid Strength = %f, Objective Function Value = %f, Computation Time (sec) = %f!\n', ...
        n, strength, cur_objective, adjoint_time);
    
    eps(3:end-2,3:end-2) = eps(3:end-2,3:end-2) - step_size * (gradient_adj(3:end-2,3:end-2));
    eps(3:end-2,3:end-2) = min(eps(3:end-2,3:end-2), max_eps_array(3:end-2,3:end-2));
    eps(3:end-2,3:end-2) = max(eps(3:end-2,3:end-2), min_eps_array(3:end-2,3:end-2));
    
end

z = reshape(apply_new_func(eps(:), strength, center), size(eps));
eps_binary = min_eps_value * (eps <= eps_midpt) + max_eps_value * (eps > eps_midpt);

% Simulate the binarized result.
simulate(spec, eps_binary, [simulation_width simulation_height]);
