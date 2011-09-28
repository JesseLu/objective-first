function [eps] = mode_cutout(eps, side)

[x, y] = border_find(side);

% Cut out the structure.
eps.x = eps.x(x,y);
eps.y = eps.y(x,y);
