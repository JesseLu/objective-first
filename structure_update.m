function [phi] = structure_update(phi, phys_res, grad_res)

    %
    % Calculate the derivative.
    %

g = real(grad_res(phi));

% Move phi with respect to the gradient.
phi = phi - g / max(abs(g(:)));

phys_res(phi)

plot_fields(size(phi), {'grad', grad_res(phi)});
lset_plot(phi)
