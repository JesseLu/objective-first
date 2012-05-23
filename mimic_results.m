function mimic_results()
    dims = [40 40];
    eps_init = 9.0;
    eps_lims = [1 12.25];

%     %% Object mimic.
%     % Set up simulation which will find what we need to mimic.
%     pad = 60;
%     ext_dims = dims + [2*pad 0];
%     ext_eps = ones(ext_dims);
%     ext_eps = my_circle(ext_eps, ext_dims/2, 8, -2);
% 
%     % Set up spec.
%     eps0 = eps_init * ones(dims);
%     eps0([1:2, end-1:end], :) = 1.0;
%     spec = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
% 
%     % Simulate and hack the spec.
%     [Ex, Ey, Hz0] = ob1_fdfd(spec.omega, ext_eps, spec.in, spec.bc);
%     figure(2)
%     ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Hz|', abs(Hz0)}, {'Re(Hz)', real(Hz0)});
%     figure(1)
%     spec.Hz0([1:2, end-1:end], :) = Hz0([pad+[1:2], end-1-pad:end-pad], :);
% 
%     % Solve.
%     eps = solve(spec, 4, 1e-6);
% 
%     % Simulate the mimic.
%     ext_eps(pad+1:pad+dims(1), :) = eps;
%     [Ex, Ey, Hz1] = ob1_fdfd(spec.omega, ext_eps, spec.in, spec.bc);
%     ob1_plot(size(ext_eps), {'\epsilon', ext_eps}, {'|Hz|', abs(Hz1)}, {'Re(Hz)', real(Hz1)});
% 
%     % Compute error.
%     comp_reg = ones(ext_dims);
%     comp_reg(pad+1:pad+dims(1), :) = 0;
%     field_err = comp_reg .* (Hz0 - Hz1);
%     err = norm(field_err) / norm(comp_reg(:) .* Hz0(:))
% 
%     ob1_plot(size(ext_eps), {'|Hz|', abs(field_err)}, {'Re(Hz)', real(field_err)});


    %% Negative index mimic.
    % Set up simulation which will find what we need to mimic.
    dims = [60 80];
    pad = 100;
    width = 4; % Of source.
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
    eps = solve(spec, 400, 1e-6);

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


