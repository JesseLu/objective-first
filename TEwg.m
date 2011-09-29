function [epsilon] = TEwg(l, w, eps, shape)


ny = floor(shape(2)/2);
epsilon = [ repmat(my_square(l, w, eps), 1, ny), ...
            my_square(l, 0, eps), ...
            repmat(my_square(l, w, eps), 1, ny)];

epsilon = repmat(epsilon, shape(1), 1);
            

function [z] = my_square(l, w, vals)

z = vals(2) * ones(l);
z(1+w:end-w, 1+w:end-w) = vals(1);
