function [epsilon] = make_epsilon(x_start, a, r, eps, dims)

% Create a waveguide mode out of a square-lattice photonic crystal.
epsilon = eps(2) * ones(dims);
hole = @(epsilon, x, y) my_hole([x, y], r, eps(1), epsilon);
x_start = x_start + a/2 + 1;
for i = 0 : dims(1)/a 
    for j = 1 : 0.5 * dims(2)/a + 1
    epsilon = hole(epsilon, x_start + i * a, (dims(2)+1)/2 + j * a);
    epsilon = hole(epsilon, x_start + i * a, (dims(2)+1)/2 - j * a);
    end
end


function [epsilon] = my_hole(pos, r, val, epsilon)

dims = size(epsilon);
[x, y] = ndgrid(1:dims(1), 1:dims(2));

z = sqrt((x - pos(1)).^2 + (y - pos(2)).^2) - r;

% Invert so that intermediary values are between 0 and 1.
z = 0.5 * (1 - z);

% Cap z.
z = 0 * (z <= 0) + 1 * (z >= 1) + z .* ((z > 0) & (z < 1));

epsilon = epsilon + (val - epsilon) .* z;
