function [phi] = structure_update(phi, phys_res, grad_res)

% alpha = -0.9:0.1:0.9;
%     g = real(grad_res(phi));
% for k = 1 : length(alpha)
%     phi_test = update_interface(phi, struct('x', 0, 'y', 0), g, 0, alpha(k));
%     r(k) = phys_res(phi_test);
% end
% plot(alpha, r)
% return

% Calculate the residual.
r = phys_res(phi)

alpha = 0.9;
for k = 1 : 1000
    % Calculate the derivative.
    g = real(grad_res(phi));

    while true
        % Update the interface using the values of the gradient.
        phi_test = update_interface(phi, struct('x', 0, 'y', 0), g, 0, alpha);
        if (phys_res(phi_test) < r(end))
            phi = phi_test;
            break
        else
            alpha = alpha / 2;
        end
    end

    r(end+1) = phys_res(phi);
    fprintf('%d: %e %e\n', k, alpha, r(end));

    if (2 * alpha > 0.9)
        alpha = 0.9;
    else
        alpha = 2 * alpha;
    end

    subplot 111; lset_plot(phi);
end
    % Move phi with respect to the gradient.
