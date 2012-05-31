function cloak_results()
    dims = [60 100]; % Dimensions of all the couplers.
    eps_init = 9.0; % Initial value of epsilon within design area.
    eps_lims = [1 12.25]; % Limited range of epsilon.
    num_iters = 400; % Number of iterations to run the design optimization for.

%     % Create the specifications for the various design problems.
%     [specs, blocks] = create_specs(dims, eps_lims, eps_init);
% 
%     for k = 1 : length(specs)
%         eps{k} = solve(specs{k}, num_iters, 1e-6, blocks{k});
%     end
% 
%     save('precomp_cloak_results.mat', 'eps', 'specs');

    results = load('precomp_cloak_results.mat', 'eps', 'specs');
    specs = results.specs;
    eps = results.eps;

% Add empty ball and channel results.
kill_pos = @(z) z .* (z < 0) + 1 * (z >= 0);
specs = [specs, {specs{2}}, {specs{5}}];
eps = [eps, {kill_pos(eps{2})}, {kill_pos(eps{5})}];
% Generate the figures.
try
    system('mkdir fig');
end
% for k = 1 : length(specs)
for k = 5:7
    basename = ['fig/cloak', num2str(k), '/'];
    try
        system(['mkdir ', basename]);
    end
    fprintf('\n\nGenerating plots used for result #%d...\n===\n', k);

    fprintf('Simulation results:\n'); % Simulate.
    [eff(k), eps_sim, Ex, Ey, Hz] = simulate(specs{k}, eps{k}, [160 100]);

    % fprintf('\nInput/output mode profiles (figure 1)...\n');
    figure(1); subplot 111; % Generate image files.
    ob1_area_plot(specs{k}.in.Hz, [basename, 'a']);
    ob1_area_plot(specs{k}.out.Hz, [basename, 'b']);
    figure(1); % Plot images for user.
    ob1_plot_images({[basename, 'a.png'], 'Input mode (Hz)'}, ...
                   {[basename, 'b.png'], 'Output mode (Hz)'});
                   
    % fprintf('\nDesign results (figure 2)...\n');
    figure(2); subplot 111; % Generate image files.
    if (min(eps{k}(:)) < 0) % Custom colormap if we have metallic devices.
        cmap = flipud([colormap('bone'); flipud(fliplr(colormap('bone')))]);


        min_eps = -2; % Override.
        max_eps = 12.25; % Override.

        % Ratio (centers around 1).
        r = (abs(min_eps) + 1) / (max_eps - 1); 
        n = size(cmap, 1);
        ind = round(n/2 * (1 - r))+2;
        cmap = interp1(cmap, ind:(n-ind)/64:n);
    else
        min_eps = min(eps{k}(:));
        max_eps = max(eps{k}(:));
        cmap = flipud(colormap('bone'));
    end
    c = 3; % Controls the color range.
    ob1_imagesc(eps_sim, cmap, ...
        [min_eps, max_eps], [basename, 'c']);
    ob1_imagesc(abs(Hz), colormap('hot'), ...
        c * mean(max(max(abs(Hz(1:20,:))))) * [0 1], [basename, 'd']);
    ob1_imagesc(real(Hz), colormap('jet'), ...
        c * mean(max(max(abs(real(Hz(1:20,:)))))) * [-1 1], [basename, 'e']);
    figure(2); % Plot images for user.
    ob1_plot_images({[basename, 'c.png'], 'Relative permittivity'}, ...
                   {[basename, 'd.png'], '|Hz| (simulation)'}, ...
                   {[basename, 'e.png'], 'Re(Hz) (simulation)'});

    fprintf('\nPress enter to continue...\n'); pause; % Wait for user.
end


function [specs, blocks] = create_specs(dims, eps_lims, eps_init)

	% Anti-reflection coating.
	eps0 = ones(dims);
	eps0(3:end,:) = 12.25;
	specs{1} = setup(0.10, eps0, eps_lims, [1 1], 'periodic');
    blocks{1} = false * ones(dims);

	% Wrap cloak. 
	block = false * ones(dims);
	block = my_circle(block, dims/2, 10, 1);
	eps0 = ones(dims);
	eps0(3:end-2,:) = eps_init;
	eps0 = my_circle(eps0, dims/2, 10, 1);
	eps0 = my_circle(eps0, dims/2, 8, -2);
	specs{2} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    blocks{2} = block;
	
	% Front/back cloak.
	block = false * ones(dims);
	block(dims(1)/2 + [-10:10],:) = true;
	eps0 = ones(dims);
	eps0 = my_circle(eps0, dims/2, 8, -2);
	eps0 = block .* eps0 + ~block .* eps_init;
	eps0([1:2, end-1:end],:) = 1;
	specs{3} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    blocks{3} = block;
	
	% Side cloak.
	block = false * ones(dims);
	block(:,dims(2)/2 + [-10:10]) = true;
	eps0 = ones(dims);
	eps0 = my_circle(eps0, dims/2, 8, -2);
	eps0 = block .* eps0 + ~block .* eps_init;
	eps0([1:2, end-1:end],:) = 1;
	specs{4} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    blocks{4} = block;

	% Channeler
	block = false * ones(dims);
	block(dims(1)/2 + [-4:4],:) = true;
	eps0 = ones(dims);
	eps0(dims(1)/2 + [-2:2], :) = -20;
	eps0(dims(1)/2 + [-2:2], dims(2)/2 + [-1:1]) = 1;
	eps0 = block .* eps0 + ~block .* eps_init;
	eps0([1:2, end-1:end],:) = 1;
	specs{5} = setup(0.15, eps0, [1 12.25], [1 1], 'periodic');
    blocks{5} = block;




function [z] = my_circle(z, pos, radius, val)
    [x, y] = ndgrid(1:size(z,1), 1:size(z,2));
    r = sqrt((x-pos(1)).^2 + (y-pos(2)).^2);
    in = r <= radius;
    z = in * val + ~in .* z;


