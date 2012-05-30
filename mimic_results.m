function mimic_results()
%     num_iters = 400;
%     [eps{1}, Hz0{1}, Hz1{1}] = object(num_iters);
%     [eps{2}, Hz0{2}, Hz1{2}] = negative_index(num_iters);
%     [eps{3}, Hz0{3}, Hz1{3}] = lens(35, 100, num_iters); % "Wide" focus.
%     [eps{4}, Hz0{4}, Hz1{4}] = lens(12, 50, num_iters); % Tight focus.
%     [eps{5}, Hz0{5}, Hz1{5}] = lens(12, 150, num_iters); % Tight focus.
%     [eps{6}, Hz0{6}, Hz1{6}] = litho(6, 3.5, 15, num_iters);
% 
%     save('precomp_mimic_results.mat', 'eps', 'Hz0', 'Hz1');

    results = load('precomp_mimic_results.mat', 'eps', 'Hz0', 'Hz1');
    eps = results.eps;
    Hz0 = results.Hz0;
    Hz1 = results.Hz1

    % Comparison lines.
    comp_lines = [300, 230, 320, 270, 370, 220];
    black_out = [220, 230, 220, 220, 220, 220];
    
% Generate the figures.
try
    system('mkdir fig');
end
for k = 1 : length(eps)
    basename = ['fig/mimic', num2str(k), '/'];
    try
        system(['mkdir ', basename]);
    end
    fprintf('\n\nGenerating plots used for result #%d...\n===\n', k);

    % Get the comparison line.
    target = Hz0{k}(comp_lines(k), :);
    actual = Hz1{k}(comp_lines(k), :);

    figure(1); subplot 111;
    ob1_area_plot(abs(target(:)), [basename, 'dt'], 'pos');
    ob1_area_plot(abs(actual(:)), [basename, 'da'], 'pos');
    figure(1)
    ob1_plot_images({[basename, 'dt.png'], 'Target (Hz)'}, ...
                   {[basename, 'da.png'], 'Actual (Hz)'});

    % Calculate the error.
    t = abs(target(:)) / norm(target(:));
    a = abs(actual(:)) / norm(actual(:));
    fprintf('Relative error on comparison line: %1.5f%%\n', 100 * norm(t - a) / norm(t));

    % Pad epsilon.
    pad = (size(Hz0{k}, 1) - size(eps{k}, 1)) / 2;
    eps_sim = [ones(pad, size(eps{k}, 2)); eps{k}; ones(pad, size(eps{k}, 2))];

    % Trim.
    trim_start = 151;
    eps_sim = eps_sim(trim_start:end,:);
    Hz0{k}(1:219,:) = 0; % Black-out field to the left of the device.
    Hz0{k} = Hz0{k}(trim_start:end, :);
    Hz1{k} = Hz1{k}(trim_start:end, :);

    c = 1.5; % Controls the color range.
    ob1_imagesc(eps_sim, colormap('bone'), ...
        [min(eps{k}(:)), max(eps{k}(:))], [basename, 'a']);

    xline = comp_lines(k) - trim_start + 1;
    ob1_imagesc(abs(Hz0{k}), colormap('hot'), ...
        c * max(abs(target(:))) * [0 1], [basename, 'b'], xline);
    ob1_imagesc(real(Hz0{k}), colormap('jet'), ...
        c * max(real(target(:))) * [-1 1], [basename, 'c'], xline);
    ob1_imagesc(abs(Hz1{k}), colormap('hot'), ...
        c * max(abs(actual(:))) * [0 1], [basename, 'e'], xline);
    ob1_imagesc(real(Hz1{k}), colormap('jet'), ...
        c * max(real(actual(:))) * [-1 1], [basename, 'f'], xline);
    figure(2)
    ob1_plot_images({[basename, 'a.png'], 'Relative permittivity'}, ...
                   {[basename, 'b.png'], '|Hz| (target)'}, ...
                   {[basename, 'c.png'], 'Re(Hz) (target)'}, ...
                   {[basename, 'e.png'], '|Hz| (simulation)'}, ...
                   {[basename, 'f.png'], 'Re(Hz) (simulation)'});

    fprintf('\nPress enter to continue...\n'); pause; % Wait for user.
end

;

function [eps, Hz0, Hz1] = object(num_iters)
    %% Object mimic.
    % Set up simulation which will find what we need to mimic.
    dims = [40 120];
    eps_init = 9.0;
    eps_lims = [1 12.25];

    ext_dims = [400 120]
    pad = 180;
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
    eps = solve(spec, num_iters, 1e-6);

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


function [eps, Hz0, Hz1] = negative_index(num_iters)
    %% Negative index mimic.
    % Set up simulation which will find what we need to mimic.
    eps_init = 9.0;
    eps_lims = [1 12.25];
    dims = [60 120];
    pad = 170;
    width = 6; % Of source.
    ext_dims = dims + [2*pad 0];

    y_pos = [1 : dims(2)] - dims(2)/2;
    gauss = @(pos, sigma) exp(-(y_pos - pos).^2 ./ (2*sigma));

    source = zeros(ext_dims);
    source(pad+dims(1)-1, :) = -gauss(0, width);
    source(pad+dims(1), :) = +gauss(0, width);

    % Simulate the negative-index material.
    [Ex, Ey, Hz0] = ob1_fdfd_adv(0.15, {ones(ext_dims), ones(ext_dims)}, ...
                                ones(ext_dims), source, 'per', 20);
    figure(1)
    ob1_plot(ext_dims, {'|Hz|', abs(Hz0)}, {'Im(Hz)', imag(Hz0)});
    % ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Ey|', abs(Ey)}, {'Re(Ey)', real(Ey)});
%     plot([  real(Hz0(pad-12+2,:)); ...
%             imag(Hz0(pad-12+2,:)); ...
%             abs(Hz0(pad-12+2,:))]', '.-');
    
    % Set up spec.
    spec.Hz0 = zeros(dims);
    spec.Hz0(1:2, :) = Hz0(pad+dims(1)+1:pad+dims(1)+2,:);
    spec.Hz0(end-1:end, :) = spec.Hz0(1:2, :);
    spec.bc = 'per';
    spec.eps0 = eps_init * ones(dims);
    spec.eps0([1:2, end-1:end],:) = 1;
    spec.eps_lims = eps_lims;
    spec.omega = 0.15;

    % Solve.
    eps = solve(spec, num_iters, 1e-6);

    % Simulate the mimic.
    ext_eps = ones(ext_dims);
    ext_eps(pad+[1:dims(1)],:) = eps;

    source = zeros(ext_dims);
    source(pad-1, :) = -gauss(0, width);
    source(pad, :) = +gauss(0, width);

    [eps_x, eps_y] = ob1_interp_eps(ext_eps); % Obtain x- and y- components of eps.
    [Ex, Ey, Hz1] = ob1_fdfd_adv(0.15, {eps_x, eps_y}, ...
                                ones(ext_dims), source, 'per', 20);
    ob1_plot(size(ext_eps), {'|Hz|', abs(Hz1)}, {'Re(Hz)', real(Hz1)});

    % Compute error.
    comp_reg = ones(ext_dims);
    comp_reg(1:pad+dims(1), :) = 0;
    field_err = comp_reg .* (Hz0 - Hz1);
    err = norm(field_err) / norm(comp_reg(:) .* Hz0(:))

    figure(1)
    ob1_plot(size(ext_eps), {'|Hz0|', abs(comp_reg .* Hz0)}, ...
                            {'|Hz1|', abs(comp_reg .* Hz1)}, ...
                            {'|err|', abs(field_err)});
    figure(2)
    ob1_plot(size(ext_eps), {'|Hz0|', abs(Hz0)}, ...
                            {'|Hz1|', abs(Hz1)}, ...
                            {'|err|', abs(field_err)});
    
function [z] = my_circle(z, pos, radius, val)
    [x, y] = ndgrid(1:size(z,1), 1:size(z,2));
    r = sqrt((x-pos(1)).^2 + (y-pos(2)).^2);
    in = r <= radius;
    z = in * val + ~in .* z;


function [eps, Hz0, Hz1] = lens(fwhm, d, num_iters)
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


function [eps, Hz0, Hz1] = litho(res, cc, d, num_iters)
    dims = [40 120];
    ext_dims = [400 120];
    eps_init = 9.0;
    eps_lims = [1 12.25];
    omega = 0.25;
    t_pml = 20;
    loc = ext_dims(1)/2 + dims(1)/2; % Location where device ends.

    % Set up simulation which will find what we need to mimic.
    y_pos = [1 : ext_dims(2)] - ext_dims(2)/2;

    target = 0;
    for k = -1 : 1
        target = target + (((y_pos + k * res) < res/4) & ((y_pos + k * res) >= -res/4));
    end
    % target = target - mean(target); % Take away DC component.

    t2  = ob1_backtrack(omega, target, 0, cc * omega);
    % t2 = t2 - mean(t2);
%     plot(([target; t2])', '.-');
%     return



    source = zeros(ext_dims);
    source(loc-d,:) = ob1_backtrack(omega, target, d, cc * omega) + 0;

    [Ex, Ey, Hz0] = ob1_fdfd_adv(omega, {ones(ext_dims), ones(ext_dims)}, ...
                                ones(ext_dims), source, 'per', t_pml, 'force');
    power = sum(conj(Ey) .* Hz0, 2); % Power calculation.

%     ob1_plot(ext_dims, {'|Hz|', abs(Hz0)}, {'Re(Hz)', real(Hz0)});
%     subplot 212; plot(abs(power), '.-');

    mean(source(loc-d,:))
    subplot 111; plot(abs([target; t2; Hz0(loc,:)]'), '.-');
    % pause


    power = abs(power(end-t_pml-1));
    Ex = power.^-0.5 * Ex;
    Ey = power.^-0.5 * Ey;
    Hz0 = power.^-0.5 * Hz0;

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

    subplot 311; plot(abs([Hz0(loc,:); Hz1(loc,:)]'), '.-');

%     figure(2);
%     ob1_plot(size(ext_eps), {'|Hz|', abs(field_err)}, {'Re(Hz)', real(field_err)});


