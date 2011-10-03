function ob1_sim(omega, epsilon, mode_num, varargin)
path('sim', path)
dims = size(epsilon);

% Expansion factor.
if ~isempty(varargin)
    c = varargin{1}; 
else
    c = [1 1];
end


% Expand epsilon.
epsilon = cat(1, repmat(epsilon(1,:), round(c(1)*dims(1)/2), 1), epsilon, ...
            repmat(epsilon(end,:), round(c(1)*dims(1)/2), 1));
epsilon = cat(2, repmat(epsilon(:,1), 1, round(c(2)*dims(2)/2)), epsilon, ...
            repmat(epsilon(:,end), 1, round(c(2)*dims(2)/2)));
dims = size(epsilon);
n = numel(epsilon);

[A, S] = ob1_matrices(dims, 0);
eps = A{2} * epsilon(:);
eps = struct('x', reshape(eps(1:n), dims), 'y', reshape(eps(n+1:2*n), dims));

[Ex, Ey, Hz] = sim_epsilon(omega, eps, '-x', mode_num);

figure(3); 
plot_fields(dims, {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});
% plot_fields(dims, {'|Ex|', real(Ex)}, {'|Ey|', real(Ey)}, {'|Hz|', real(Hz)});

% Dx = eps.x .* Ex;
% Dy = eps.y .* Ey;
% figure(4);
% plot_fields(dims, {'|Dx|', abs(Dx)}, {'|Dy|', abs(Dy)}, {'|Hz|', abs(Hz)});


% Extra plotting for postdeadline Fios submission.
figure(4);
dims = size(epsilon)
plot_fields([160 100], {'', real(Ey(21:180,31:130))})

figure(5);
imagesc(reshape(real(epsilon(21:180,31:130)), [160 100])'); axis equal tight;
set(gca, 'YDir', 'normal');
cmap = colormap('bone');
colormap(cmap(end:-1:1,:));
