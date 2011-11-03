function [t] = backtrack_linesearch(f, t, f0, g0dx, alpha, beta);
% Backtracking line search on one-dimensional function f.
% Assumes alpha in (0, 0.5) and beta in (0, 1).

% y = []; % For debugging purposes.
% x = [];
while f(t) > f0 + alpha * t * g0dx
%     y(end+1) = f(t);
%     x(end+1) = t;
    t = beta * t;
    if (t <= 1e-6) % Just try to get some improvement.
        warning('Setting alpha to 0');
        alpha = 0;
    end
    if (t <= eps) % Not a descent direction.
%         semilogx(x, y - f0, '.-');
%         drawnow
        error('Backtracking line-search failed.');
        break
    end
end

