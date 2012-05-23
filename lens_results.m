function mimic_results(fwhm, d)
    dims = [40 120];
    ext_dims = [400 120];
    eps_init = 9.0;
    eps_lims = [1 12.25];
    omega = 0.25;
    t_pml = 20;

    % Set up simulation which will find what we need to mimic.
    y_pos = [1 : ext_dims(2)] - ext_dims(2)/2;
    gauss = @(pos, sigma) exp(-(y_pos - pos).^2 ./ (2*sigma^2));
    h = ob1_backtrack(omega, gauss(0, fwhm / (2*sqrt(2*log(2)))), 0, omega);

    source = zeros(ext_dims);
    source(t_pml+1,:) = -i * omega * ob1_backtrack(omega, h, ext_dims(1)/2 + dims(1)/2 + d - t_pml, omega);

    [Ex, Ey, Hz0] = ob1_fdfd_adv(omega, {ones(ext_dims), ones(ext_dims)}, ...
                                ones(ext_dims), source, 'per', t_pml);

    % Set up spec.
    


    ob1_plot(ext_dims, {'|Hz|', abs(Hz0)}, {'Re(Hz)', real(Hz0)});
    return
    %% Object mimic.
    pad = 60;
    ext_dims = dims + [2*pad 0];
    ext_eps = ones(ext_dims);
    ext_eps = my_circle(ext_eps, ext_dims/2, 8, -2);

    % Set up spec.
    eps0 = eps_init * ones(dims);
    eps0([1:2, end-1:end], :) = 1.0;
    spec = setup(0.15, eps0, eps_lims, [1 1], 'periodic');

    % Simulate and hack the spec.
    [Ex, Ey, Hz0] = ob1_fdfd(spec.omega, ext_eps, spec.in, spec.bc);
    figure(2)
    ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Hz|', abs(Hz0)}, {'Re(Hz)', real(Hz0)});
    figure(1)
    spec.Hz0([1:2, end-1:end], :) = Hz0([pad+[1:2], end-1-pad:end-pad], :);

    % Solve.
    eps = solve(spec, 4, 1e-6);

    % Simulate the mimic.
    ext_eps(pad+1:pad+dims(1), :) = eps;
    [Ex, Ey, Hz1] = ob1_fdfd(spec.omega, ext_eps, spec.in, spec.bc);
    ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Hz|', abs(Hz1)}, {'Re(Hz)', real(Hz1)});

    % Compute error.
    comp_reg = ones(ext_dims);
    comp_reg(pad+1:pad+dims(1), :) = 0;
    field_err = comp_reg .* (Hz0 - Hz1);
    err = norm(field_err) / norm(comp_reg(:) .* Hz0(:))

    ob1_plot(size(ext_eps), {'|Hz|', abs(field_err)}, {'Re(Hz)', real(field_err)});


