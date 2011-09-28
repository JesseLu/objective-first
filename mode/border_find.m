function [x, y] = border_find(side)

% global MY_DEVICE
% device = MY_DEVICE;
% if ~isempty(device)
%     % Device boundaries.
%     lim.x = device.offset(1) + [0 device.dims(1)-1];
%     lim.y = device.offset(2) + [0 device.dims(2)-1];
% 
% else
    % Get device location and size.
    global DIMS_ 
    dims = DIMS_;

    % Device boundaries.
    lim.x = [1 dims(1)];
    lim.y = [1 dims(2)];
% end

    %
    % Find the relevant slice of the structure.
    %

% Determine where to "cut out" the structure.
switch (side)
    case 'x-'
        x = lim.x(1);
        y = lim.y(1) : lim.y(2);
    case 'x+'
        x = lim.x(2);
        y = lim.y(1) : lim.y(2);
    case 'y-'
        x = lim.x(1) : lim.x(2);
        y = lim.y(1);
    case 'y+'
        x = lim.x(1) : lim.x(2);
        y = lim.y(2);
end


