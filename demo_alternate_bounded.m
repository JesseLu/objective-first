function demo_alternate_bounded(dims, cgo_iters, savefile)

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
% phi2 = lset_union(phi, (eta)); % filled.
phi2 = lset_intersect(phi, lset_complement(eta));
phi = lset_complement(phi);
eta2 = lset_box([0 0], dims/2 + 2);
small_box = lset_box([0 0], dims/2 - 2);
phi2 = lset_union(phi2, lset_intersect(small_box, lset_checkered));
phi2 = lset_complement(phi2);

% Initialize phi, and create conversion functions.
[phi2p, phi2eps, phi_smooth] = setup_levelset(phi, 1.0, 12.25, 1e-3);


    %
    % Setup the constrained gradient optimization.
    %

% Objective function and its gradient.
[f, g] = em_physics(omega); 


% This constraint function allows both variables to change.
bound = @(x) 1.0 * (x < 1.0) + 12.25 * (x > 12.25) + ...
    x .* ((x >= 1.0) & (x <= 12.25));
c = @(v, dv, s) struct( 'x', v.x - s * (field_template .* dv.x), ...
                        'p', real(bound(v.p - s * ((eta2(:) < 0) .* dv.p))));

% Initial values.
[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));
v.x = [Ex(:); Ey(:); Hz(:)];
% randn('state', 1);
% v.x = randn(size(v.x));
v.p = phi2p(phi2);
v.p = v.p(:);
% v.p = randn(N, 1);


    %
    % Optimize using the c-go package.
    %

tic;
fval = [];
ss_hist = {[], []};
for k = 1 : ceil(cgo_iters/1e2)
    [v, fval0, ss_hist0] = cgo_opt(f, g, c, v, 1e2, 2.^[-20:20]); 
    fval = [fval, fval0];
    ss_hist{1} = [ss_hist{1}, ss_hist0{1}];
    ss_hist{2} = [ss_hist{2}, ss_hist0{2}];
    my_plot(v, fval, ss_hist);
    save(savefile, 'v', 'fval', 'ss_hist');
    fprintf('%d: %e\n', k*1e2, fval(end));
end
fprintf('%e, ', f(v)); toc

figure(2); quick_plot(v, fval, ss_hist{1});
figure(3); cgo_visualize(fval, ss_hist{2});

function my_plot(v, fval, ss_hist)
    
global DIMS_
dims = DIMS_;
N = prod(dims);
% figure(1); 
plot_fields(dims, ...
    {'Re(Ey)', real(v.x(N+1:2*N))}, {'|Ey|', abs(v.x(N+1:2*N))}, {'p', v.p});

% figure(2); cgo_visualize(fval, ss_hist);

drawnow



