function waves(omega, d, dims, fwhm)
    % d = 150;
    cutoff = 1 * omega;
    % width = 20;

    y_pos = [1 : dims(2)] - dims(2)/2;
    gauss = @(pos, sigma) exp(-(y_pos - pos).^2 ./ (2*sigma^2));
    h = ob1_backtrack(omega, gauss(0, fwhm / (2*sqrt(2*log(2)))), 0, cutoff);

    source = zeros(dims);
    source(170-d,:) = -i * omega * ob1_backtrack(omega, h, d, cutoff);
    % source(70-(d-1),:) = backtrack(omega, h, d-1, omega);

    [Ex, Ey, Hz0] = ob1_fdfd_adv(omega, {ones(dims), ones(dims)}, ...
                                ones(dims), source, 'per', 20);
    ob1_plot(dims, {'|Hz|', abs(Hz0)}, {'real(Hz)', real(Hz0)});

    % Power calculation.
    power = sum(conj(Ey) .* Hz0, 2);
    subplot 212; plot(abs(power), '.-');

    % Try to calculate Ey.
    Ey1 = -i/omega * diff(Hz0, 1, 1);
    
    ob1_plot(dims - [1 0], {'|Ey|', abs(Ey(2:end,:))}, {'|Ey1|', abs(Ey1)});




