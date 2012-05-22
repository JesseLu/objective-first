function mimic_results()
    dims = [40 40];
    eps_init = 9.0;
    eps_lims = [1 12.25];

    % Object mimic.
    % Set up simulation which will find what we need to mimic.
    pad = 60;
    ext_dims = dims + [2*pad 0];
    ext_eps = ones(ext_dims);
    ext_eps = my_circle(ext_eps, ext_dims/2, 8, -2);

    % Set up spec.
    eps0 = eps_init * ones(dims);
    eps0([1:2, end-1:end], :) = 1.0;
    spec = setup(0.15, eps0, eps_lims, [1 1], 'periodic');

    % Simulate and hack the spec.
    [Ex, Ey, Hz] = ob1_fdfd(spec.omega, ext_eps, spec.in, spec.bc);
    figure(2)
    ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});
    figure(1)
    spec.Hz0([1:2, end-1:end], :) = Hz([pad+[1:2], end-1-pad:end-pad], :);

    % Solve.
    eps = solve(spec, 10, 1e-6);

    [eff, eps_sim, Ex, Ey, Hz] = simulate(spec, eps, [160 100]);

    
function [z] = my_circle(z, pos, radius, val)
    [x, y] = ndgrid(1:size(z,1), 1:size(z,2));
    r = sqrt((x-pos(1)).^2 + (y-pos(2)).^2);
    in = r <= radius;
    z = in * val + ~in .* z;

