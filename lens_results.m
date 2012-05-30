function mimic_results(fwhm, d, num_iters)
    dims = [40 120];
    ext_dims = [400 120];
    eps_init = 9.0;
    eps_lims = [1 12.25];
    omega = 0.25;
    t_pml = 20;
    loc = ext_dims(1)/2 + dims(1)/2; % Location where device ends.

    % Set up simulation which will find what we need to mimic.
    y_pos = [1 : ext_dims(2)] - ext_dims(2)/2;
    gauss = @(pos, sigma) exp(-(y_pos - pos).^2 ./ (2*sigma^2));
    h = ob1_backtrack(omega, gauss(0, fwhm / (2*sqrt(2*log(2)))), 0, omega);

    source = zeros(ext_dims);
    source(ext_dims(1)/2,:) = -i * omega * ob1_backtrack(omega, h, dims(1)/2 + d, omega);

    [Ex, Ey, Hz0] = ob1_fdfd_adv(omega, {ones(ext_dims), ones(ext_dims)}, ...
                                ones(ext_dims), source, 'per', t_pml);
    power = sum(conj(Ey) .* Hz0, 2); % Power calculation.
    power = abs(power(loc));
    Ex = power.^-0.5 * Ex;
    Ey = power.^-0.5 * Ey;
    Hz0 = power.^-0.5 * Hz0;

    ob1_plot(ext_dims, {'|Hz|', abs(Hz0)}, {'Re(Hz)', real(Hz0)});
    subplot 122; plot(abs(power), '.-');
%     pause

    % Set up spec.
    eps0 = eps_init * ones(dims);
    eps0([1:2, end-1:end], :) = 1.0;
    spec = setup(omega, eps0, eps_lims, [1 1], 'periodic');
    spec.Hz0([end-1:end], :) = Hz0(loc + [-1 0], :);
    close(gcf);

    % Solve.
    eps = solve(spec, num_iters, 1e-6);

    % Simulate the mimic.
    ext_eps = ones(ext_dims);
    ext_eps(ext_dims(1)/2-dims(1)/2+[1:dims(1)], :) = eps;
    [Ex, Ey, Hz1] = ob1_fdfd(spec.omega, ext_eps, spec.in, spec.bc);
    ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Hz|', abs(Hz1)}, {'Re(Hz)', real(Hz1)});

    % Compute error.
    comp_reg = ones(ext_dims);
    comp_reg(1:loc, :) = 0;
    field_err = comp_reg .* (abs(Hz0) - abs(Hz1));
    err = norm(field_err) / norm(comp_reg(:) .* abs(Hz0(:)))

    subplot 311; plot(abs([Hz0(loc+d,:); Hz1(loc+d,:)]'), '.-');

    % ob1_plot(size(ext_eps), {'|Hz|', abs(field_err)}, {'Re(Hz)', real(field_err)});


