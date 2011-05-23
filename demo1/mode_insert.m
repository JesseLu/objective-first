function [Ex, Ey, Hz] = mode_insert(mode, side)

[x, y] = border_find(side);

global MY_DIMS

Ex = zeros(MY_DIMS);
Ey = zeros(MY_DIMS);
Hz = zeros(MY_DIMS);

switch (side)
    case '-x'


end
