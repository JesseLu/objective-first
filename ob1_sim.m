function ob1_sim(omega, epsilon)

dims = size(epsilon);

[A, S] = ob1_matrices(dims, 0);

% Expand epsilon.
epsilon = cat(1, repmat(epsilon(1,:), dims(1)/2, 1), epsilon, ...
            repmat(epsilon(end,:), dims(1)/2, 1));
epsilon = cat(2, repmat(epsilon(:,1), 1, dims(2)/2), epsilon, ...
            repmat(epsilon(:,end), 1, dims(2)/2));

imagesc(epsilon')
