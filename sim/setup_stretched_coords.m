function [A] = stretched_coords(dims, end_pos, shift, dir, thickness, ...
    sigma, order)


    %
    % Find the position of each point relative to the PML interface.
    %

switch dir
    case 'x'
        pos = layer_num([1:dims(1)]+shift(1), end_pos, thickness);
        A = repmat(pos', 1, dims(2));
    case 'y'
        pos = layer_num([1:dims(2)]+shift(2), end_pos, thickness);
        A = repmat(pos, dims(1), 1);
    otherwise
        error('Invalid DIR.');
end


    %
    % Determine sigma values.
    %

A = (ones(dims) + i * sigma * (A/thickness).^order).^-1;
A = spdiags(A(:), 0, prod(dims), prod(dims));



function [pos] = layer_num(z, end_pos, thickness)
% Calculate the position of each layer relative to the beginning of the PML.
pos = thickness - min(cat(3, [z-end_pos(1)], [end_pos(2)-z]), [], 3);
pos = pos .* (pos > 0);

