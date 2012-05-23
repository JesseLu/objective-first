function mimic_results(fwhm, f_len, num_iters)
    dims = [40 120];
    eps_init = 9.0;
    eps_lims = [1 12.25];
    omega = 0.25

    %% Object mimic.

    % Set up target.
    y_pos = [1 : dims(2)] - dims(2)/2;
    gauss = @(pos, sigma) exp(-(y_pos - pos).^2 ./ (2*sigma^2));
    target = ob1_backtrack(omega, gauss(0, fwhm / (2*sqrt(2*log(2)))), 0, omega);

    % Set up spec.
    eps0 = eps_init * ones(dims);
    eps0([1:2, end-1:end], :) = 1.0;
    spec = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    Hz0_out(1,:) = ob1_backtrack(omega, target, f_len+1, omega);
    Hz0_out(2,:) = ob1_backtrack(omega, target, f_len, omega);
    Ey_out = -i / omega * (Hz0_out(2,:) - Hz0_out(1,:)); % Assumes epsilon = 1 everywhere.
    Hz_out = mean(Hz0_out, 1);
    power = abs(Ey_out(:)' * Hz_out(:));
    spec.Hz0(end-1:end,:) = power^-0.5 * Hz0_out;
    Ey_out = -i / omega * (spec.Hz0(end-1,:) - spec.Hz0(end,:)); % Assumes epsilon = 1 everywhere.
    Hz_out = spec.Hz0(end,:);
    power = abs(Ey_out(:)' * Hz_out(:))
    target = power^-0.5 * target;
    power
    close(gcf);

    spec.Hz0(1,:) = ob1_backtrack(omega, target, f_len+dims(1), omega);
    spec.Hz0(2,:) = ob1_backtrack(omega, target, f_len+dims(1)-1, omega);
    % spec.eps0 = ones(dims);

    % Solve.
    eps = solve(spec, num_iters, 1e-6);

    % Simulate.
    [eff, eps_sim, Ex, Ey, Hz] = simulate(spec, eps, dims + [4*f_len 0]);

    % Calculate error
    result = Hz(dims(1) + 3*f_len,:);
    result = result(:);
    target = target(:);
    err = norm(abs(target) - abs(result)) / norm(target)


    figure(2)
    plot([abs(target), abs(result)], '.-');
