function ob1_sim(omega, epsilon)

dims = size(epsilon);


% Expand epsilon.
epsilon = cat(1, repmat(epsilon(1,:), dims(1)/2, 1), epsilon, ...
            repmat(epsilon(end,:), dims(1)/2, 1));
epsilon = cat(2, repmat(epsilon(:,1), 1, dims(2)/2), epsilon, ...
            repmat(epsilon(:,end), 1, dims(2)/2));
dims = size(epsilon);
n = numel(epsilon);

[A, S] = ob1_matrices(dims, 0);
eps = A{2} * epsilon(:);
eps = struct('x', reshape(eps(1:n), dims), 'y', reshape(eps(n+1:2*n), dims));

path('~/wave-tools/em_2dte', path);
sim_epsilon(omega, eps, '-x');
