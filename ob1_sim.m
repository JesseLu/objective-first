function ob1_sim(omega, epsilon, mode_num, varargin)
path('mode', path)
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

[Ex, Ey, Hz] = sim_epsilon(omega, eps, 'x-', mode_num);

figure(3); 
plot_fields(dims, {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});
plot_fields(dims, {'|Ey|', abs(Ey)});
% plot_fields(dims, {'|Ex|', real(Ex)}, {'|Ey|', real(Ey)}, {'|Hz|', real(Hz)});

% Dx = eps.x .* Ex;
% Dy = eps.y .* Ey;
% figure(4);
% plot_fields(dims, {'|Dx|', abs(Dx)}, {'|Dy|', abs(Dy)}, {'|Hz|', abs(Hz)});

