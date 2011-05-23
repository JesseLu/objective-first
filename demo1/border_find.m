function [x, y] = border_find(side)

% Get device location and size.
global MY_DEVICE
device = MY_DEVICE;


    %
    % Find the relevant slice of the structure.
    %

% Device boundaries.
lim.x = device.offset(1) + [0 device.dims(1)-1];
lim.y = device.offset(2) + [0 device.dims(2)-1];

% Determine where to "cut out" the structure.
switch (side)
    case '-x'
        x = lim.x(1) - 1;
        y = lim.y(1) : lim.y(2);
    case '+x'
        x = lim.x(2) + 1;
        y = lim.y(1) : lim.y(2);
    case '-y'
        x = lim.x(1) : lim.x(2);
        y = lim.y(1) - 1;
    case '+y'
        x = lim.x(1) : lim.x(2);
        y = lim.y(2) + 1;
end


