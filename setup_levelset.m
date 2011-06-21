function [phi2p, phi2e, smooth_phi] = ...
    setup_levelset (phi, eps_lo, eps_hi, alpha_smooth)


global DIMS_ S_ 
dims = DIMS_;
% Epsilon is defined on Ez, interpolate to get epsilon at Ex and Ey.
A_spread_x = 0.5 * [S_(0,0)+S_(1,0)]; 
A_spread_y = 0.5 * [S_(0,0)+S_(0,1)];
p2eps = @(p) struct('x', reshape(A_spread_x * p(:), dims), ...
                    'y', reshape(A_spread_y * p(:), dims));
                    % 'z', reshape(p(:), dims)); % For 3D code.
A_spread = [A_spread_x; A_spread_y];
p2e = @(p) A_spread * p(:);

% Average Ex and Ey, to get back p (Ez).
A_gather = 1/4 * [S_(0,0)+S_(-1,0), S_(0,0)+S_(0,-1)];
eps2p = @(eps) reshape(A_gather * [eps.x(:); eps.y(:)], dims);
e2p = @(e) reshape(A_gather * e, dims);

% Helper function to convert from phi to p ("digitize phi").
% Note that this is a "soft" digitization because intermediate values of
% epsilon are allowed near and on the material interfaces.
% phi2p = @(phi) ((phi .* (phi > -1 & phi < 1) + (phi >= 1) - (phi <= -1)) * ... 
%     (eps_hi-eps_lo)/2) + (eps_hi+eps_lo)/2;
phi2p = @(phi)  lset_phi2p(phi)/2 * (eps_hi-eps_lo) + (eps_hi+eps_lo)/2;

% Make it easier to go straight to epsilon from phi.
phi2e = @(phi) p2e(phi2p(phi));
phi2eps = @(phi) p2eps(phi2p(phi));

% Re-initialization function to keep phi "well-behaved".
smooth_phi = @(phi) signed_distance(phi, alpha_smooth);

% Re-initialize phi.
phi = smooth_phi(phi);
