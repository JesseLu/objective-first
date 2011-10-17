function simple_test(x0, t)

n = length(x0);
b = randn(n, 1);
b = -1 * (b < 0) + -1 * (b >= 0);
A = randn(n);

cvx_begin
    variable x_star(n)
    minimize norm(A * x_star - b, 2)
    subject to 
        x_star >= 0
        x_star <= 1
cvx_end

types = {'ncg', 'lbfgs', 'tn'};
for k = 1 : length(types)
    hist(:,k) = my_test(A, b, x0, x_star, t, types{k});
end
subplot 111; semilogy(hist, '.-'); ylabel('mean of x');
legend(types)

function [hist] = my_test(A, b, x0, x_star, t, type)
for k = 1 : length(t)
    fprintf('\n t = %e\n', t(k));
    switch type
        case 'ncg'
            out = ncg(@(x) pob(A, b, x, t(k)), x0);
        case 'lbfgs'
            out = lbfgs(@(x) pob(A, b, x, t(k)), x0);
        case 'tn'
            out = tn(@(x) pob(A, b, x, t(k)), x0);
    end
    x0 = out.X;
    hist(k) = norm(x0 - x_star);
end
