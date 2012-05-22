function metamaterial_results()
    dims = [40 80]; % Dimensions of all the couplers.
    eps_init = 9.0; % Initial value of epsilon within design area.
    eps_lims = [1 12.25]; % Limited range of epsilon.
    num_iters = 40; % Number of iterations to run the design optimization for.

    % Create the specifications for the various design problems.
    [specs, blocks] = create_specs(dims, eps_lims, eps_init);

    for k = 1 : length(specs)
        eps{k} = solve(specs{k}, num_iters, 1e-6, blocks{k});
    end

    save('meta_results.mat', 'eps', 'specs');

    results = load('meta_results.mat', 'eps', 'specs');
    specs = results.specs;
    eps = results.eps;

    % Generate the figures.
    for k = 1 : length(specs)
        [eff(k), eps_sim, Ex, Ey, Hz] = simulate(specs{k}, eps{k}, [160 100]);
        pause
    end

function [specs, blocks] = create_specs(dims, eps_lims, eps_init)

	% Anti-reflection coating.
	eps0 = ones(dims);
	eps0(3:end,:) = 12.25;
	specs{1} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
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

