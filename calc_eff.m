function eff_calc(omega, epsilon, mode_num, varargin)
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

% Create the reference epsilon.
epsilon0 = repmat(epsilon(1,:), size(epsilon, 1), 1);

% Reference simulation.
[A, S] = ob1_matrices(dims, 0);
eps = A{2} * epsilon0(:);
eps = struct('x', reshape(eps(1:n), dims), 'y', reshape(eps(n+1:2*n), dims));
[Ex, Ey, Hz] = sim_epsilon(omega, eps, '-x', mode_num(1));
[P0mode, P0tot] = ...
    my_power(omega, eps, Ey, Hz, size(epsilon, 2) - 25, mode_num(1));

% Device simulation.
[A, S] = ob1_matrices(dims, 0);
eps = A{2} * epsilon(:);
eps = struct('x', reshape(eps(1:n), dims), 'y', reshape(eps(n+1:2*n), dims));
[Ex, Ey, Hz] = sim_epsilon(omega, eps, '-x', mode_num(1));
[Pmode, Ptot] = ...
    my_power(omega, eps, Ey, Hz, size(epsilon, 2) - 25, mode_num(2));

% Print out results.
fprintf('Ref: %1.3f out of %1.3f = %1.3f%% of power in target mode\n', ...
    P0mode, P0tot, P0mode/P0tot*100);
fprintf('Device: %1.3f out of %1.3f = %1.3f%% of power in target mode\n', ...
    Pmode, Ptot, Pmode/Ptot*100);
fprintf('Efficiency: %1.3f%%\n', Pmode/P0mode*100);




function [Pmode, Ptot] = my_power(omega, eps, E, H, x_pos, mode_num)

% Set the device boundary.
dev.dims = round(size(E) - [50 0]);
dev.offset = [25 1];
global MY_DEVICE
MY_DEVICE = dev;

% Isolate relevant slice of fields.
E = E(x_pos,:);
H = H(x_pos,:);

% Find total power.
Ptot = abs(sum(E .* H));
Ptot = abs(E(:)' * H(:));

% Find the relevant mode.
input = mode_solve(mode_cutout(eps, '+x'), omega, '-x', mode_num);

% Isolate relevant mode.
E1 = input.Et * (E(:)' * input.Et) / norm(input.Et)^2;
H1 = input.Ht * (H(:)' * input.Ht) / norm(input.Ht)^2;
% subplot 121; plot(abs([E(:), E1(:)]), '.-'); hold on; plot(abs(input.Et)); 
% subplot 122; plot(abs([H(:), H1(:)]), '.-'); hold on; plot(abs(input.Ht)); 
% hold off; pause;

% Calculate outgoing power in the mode.
Pmode = abs(sum(E1 .* H1));
Pmode = abs(E1(:)' * H1(:));




