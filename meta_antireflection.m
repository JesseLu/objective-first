function meta_antireflection()
eps_init = 9.0
% % Anti-reflection coating.
% eps0 = ones(dims);
% eps0(3:end,:) = 12.25;
% spec = setup(0.15, eps0, [1 12.25], [1 1], 'periodic');
% close(gcf);
% eps = solve(spec, 40, 1e-6);
% simulate(spec, eps, [160 80]);

% % Wrap cloak. 
% 
% dims = [40 80];
% obj_dims = [10 10];
% 
% block = false * ones(dims);
% block = my_circle(block, dims/2, 10, 1);
% 
% eps0 = ones(dims);
% eps0(3:end-2,:) = eps_init;
% eps0 = my_circle(eps0, dims/2, 10, 1);
% eps0 = my_circle(eps0, dims/2, 8, -2);
% 
% spec = setup(0.15, eps0, [1 12.25], [1 1], 'periodic');
% close(gcf);
% 
% eps = solve(spec, 40, 1e-6, block);
% 
% simulate(spec, eps, [160 80]);
% return

% Back cloak.

dims = [80 80];
obj_dims = [10 10];

block = false * ones(dims);
block(1:40,:) = true;

eps0 = ones(dims);

% eps0(25 + [1:obj_dims(1)], ...
%     (dims(2)-obj_dims(2))/2:(dims(2)+obj_dims(2))/2) = -10;
eps0 = my_circle(eps0, [30 dims(2)/2], 8, -2);
eps0 = block .* eps0 + ~block .* eps_init;
eps0(end,:) = 1;

spec = setup(0.15, eps0, [1 12.25], [1 1], 'periodic');
close(gcf);

eps = solve(spec, 400, 1e-6, block);
simulate(spec, eps, [160 80]);


function [z] = my_circle(z, pos, radius, val)
    [x, y] = ndgrid(1:size(z,1), 1:size(z,2));
    r = sqrt((x-pos(1)).^2 + (y-pos(2)).^2);
    in = r <= radius;
    z = in * val + ~in .* z;

