function [f, g] = pob(A, b, x, bar)

if nargin == 3
    f = 0.5 * norm(A * x - b)^2;
    g = A' * (A * x - b);
else
    if (any(x <= 0) | any(x >= 1))
        f = Inf;
    else
        f = 0.5 * norm(A * x - b)^2 - ...
            sum(bar * (log(x) + log(1 - x)));
    end

    g = A' * (A * x - b) - ...
        bar * (x.^-1 - (1 - x).^-1);
end

