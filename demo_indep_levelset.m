function demo_indep_levelset(dims, cgo_iters)
% Try:
% changepath % need wave-tools and level-set packages from github.com/JesseLu
% demo_indep_levelset([80 80], 20);


path(path, '~/c-go'); % Make sure we have access to c-go.
path(path, '~/level-set'); % Make sure we have access to level-set.

omega = 0.15; % Angular frequency of desired mode.


    %
    % Helper function for determining derivative matrices.
    % Also, helper global variables for prettier argument passing.
    %

global S_ D_ DIMS_ 

% Shortcut to form a derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
D_ = @(x) spdiags(x(:), 0, numel(x), numel(x));

DIMS_ = dims;
N = prod(dims);


    %
    % Make the initial structure.
    %

lset_grid(dims);
phi = lset_box([0 0], [1000 10]);
eta = lset_box([0 0], dims/2);
small_box = lset_box([0 0], dims/2 - 10);
phi2 = lset_intersect(phi, lset_complement(small_box));
phi2 = lset_complement(phi2);
phi = lset_complement(phi);
eta2 = lset_box([0 0], dims/2 + 2);
% phi2 = lset_union(phi2, lset_intersect(small_box, lset_checkered));
% phi2 = lset_complement(phi2);

% Initialize phi, and create conversion functions.
[phi2p, phi2eps, phi_smooth] = setup_levelset(phi, 1.0, 12.25, 1e-3);


    %
    % Setup the constrained gradient optimization.
    %


% % This constraint function allows both variables to change.
% c = @(v, dv, s) struct( 'x', v.x - s * (tp .* dv.x), ...
%                         'p', v.p - s * ((eta(:) < 0) .* dv.p));

c = @(phi, dphi, s) levelset_step(phi, (eta2 < 0) .* reshape(real(dphi), dims), s); 

% Initial values.
[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));
x = [Ex(:); Ey(:); Hz(:)];
p = phi2p(phi);
p = p(:);

% Get the field.
[A, b, reinsert] = em_physics_direct('field', omega, field_template, x);
x = reinsert(A(p) \ b(p));

% Back-out field components.
Ex = reshape(x(1:N), dims);
Ey = reshape(x(N+1:2*N), dims);
Hz = reshape(x(2*N+1:end), dims);

figure(1); plot_fields(dims, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});

% Objective function and its gradient.
[f, g] = em_physics_levelset_only(omega, phi2p, x); 

 
    %
    % Optimize using the c-go package.
    %

figure(2)
tic; fprintf('Gradient solve: ');
phi = phi2;
fval = [];
ss_hist = [];
for k = 1 : ceil(cgo_iters/1)
    my_plot(phi, phi2p(phi), fval, ss_hist);
    [phi, fval0, ss_hist0] = opt(f, g, c, phi, 1, 2.^[-10:20]); 
    fval = [fval, fval0];
    ss_hist = [ss_hist, ss_hist0];
end
fprintf('%e, ', f(phi)); toc
figure(3); semilogy(fval(2:2:end), '.-');
ylabel('physics residual');
xlabel('iterations');


function my_plot(phi, p, fval, ss_hist)
    
global DIMS_
dims = DIMS_;
N = prod(dims);
% figure(1); plot_fields(dims, {'p', p}, {'p', p});
% plot_fields(dims, {'p', p}, {'p', p});
subplot 121; lset_plot(phi);
drawnow

% figure(2); cgo_visualize(fval, ss_hist);




