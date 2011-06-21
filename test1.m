% TEST1
%
% Solve for the field using generic linear algebra factor-solve tool.
help test1

dims = [80 80]; % Size of the grid.
omega = 0.15; % Angular frequency of desired mode.
path(path, '~/level-set'); % Make sure we have access to level-set.
path(path, '~/wave-tools/em_bval_2dte');
path(path, '~/wave-tools/helper');

lset_grid(dims);
phi = lset_box([0 0], [1000 10]);
eta = lset_box([0 0], [40 40]);
phi = lset_complement(lset_intersect(phi, lset_complement(eta)));
lset_plot(phi);

global S_ D_ DIMS_ 

% Shortcut to form a derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
D_ = @(x) spdiags(x(:), 0, numel(x), numel(x));

DIMS_ = dims;

N = prod(dims);


% Initialize phi, and create conversion functions.
[phi2e, phi2eps, phi_smooth] = setup_levelset(phi, 1.0, 12.25, 1e-3);

[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));

plot_fields(dims, ...
    {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, {'Re(Hz)', real(Hz)}, ...
    {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});

% make the physics matrices.
% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))];  
Hcurl = [   (S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; 

A = @(e) [Ecurl, -i*omega*speye(N); i*omega*D_(e), Hcurl];
