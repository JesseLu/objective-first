function [phi, phi2p, p2e, e2p] = setup_levelset (phi, device, eps_lo, eps_hi)

global S D

% Helper function to convert from phi to p ("digitize phi").
% Note that this is a "soft" digitization because intermediate values of
% epsilon are allowed near and on the material interfaces.
phi2p = @(phi) ((phi .* (phi > -1 & phi < 1) + (phi >= 1) - (phi <= -1)) * ... 
    (eps_hi-eps_lo)/2) + (eps_hi+eps_lo)/2;

% Epsilon is defined on Ez, interpolate to get epsilon at Ex and Ey.
A_spread = 0.5 * [S(0,0)+S(1,0); S(0,0)+S(0,1)];
p2e = @(p) A_spread * p(:);

% Average Ex and Ey, to get back p (Ez).
A_gather = 1/4 * [S(0,0)+S(-1,0), S(0,0)+S(0,-1)];
e2p = @(eps) reshape(A_gather * eps(:), dims);


