function [Ex, Ey, Hz] = mode_insert(mode, side)

[x, y] = border_find(side);

global MY_DIMS

Ex = zeros(MY_DIMS);
Ey = zeros(MY_DIMS);
Hz = zeros(MY_DIMS);

switch (side)
    case '-x'
        Ey(x,y) = mode.Jt;
        Ey(x-1,y) = -mode.Jt * exp(i * mode.beta);

    case '+x'
        Ey(x,y) = mode.Jt;
        Ey(x+1,y) = -mode.Jt * exp(i * mode.beta);

    case '-y'
        Ex(x,y) = mode.Jt;
        Ex(x,y-1) = -mode.Jt * exp(i * mode.beta);

    case '+y'
        Ex(x,y) = mode.Jt;
        Ex(x,y+1) = -mode.Jt * exp(i * mode.beta);
end
