function [x] = ob1_wgmode(omega, eps, edge, in_out)
% 
% Description
%     Solve for a particular waveguide mode at the edge of the grid.

path('~/wave-tools/em_bval_2dte', path);
global DIMS_
DIMS_ = size(eps.x);
dims = size(eps.x);

% Find mode.
mode = mode_solve(mode_cutout(eps, edge), omega, edge);

% Determine phase.
switch edge(1)
    case 'x'
        phase = (mode.beta*(dims(1)-1))/2;
    case 'y'
        phase = (mode.beta*(dims(2)-1))/2;
end

switch in_out
    case 'in'
        phase_pol = -1;
    case 'out'
        phase_pol = 1;
end

% Insert single mode into empty 2D array.
[Ex, Ey, Hz] = mode_insert(mode, edge, in_out, phase_pol * phase);

% Double up on the magnetic fields.
x = Hz(:);

% Find the required phase relations.
switch edge
    case 'x-'
        shift = 1;
        dp = -phase_pol * mode.beta;
    case 'x+'
        shift = -1;
        dp = -phase_pol * mode.beta;
    case 'y-'
        shift = dims(1);
        dp = -phase_pol * mode.beta;
    case 'y+'
        shift = -dims(1);
        dp = -phase_pol * mode.beta;
end

% Implement the added phase.
ind = find(x);
x(ind+shift) = x(ind) * exp(i * dp);

clear global DIMS_

% Plot the mode.
subplot 131; plot(abs(mode.El));
subplot 132; plot(abs(mode.Et));
subplot 133; plot(abs(mode.Ht));
pause;
