function [Ex, Ey, Hz] = setup_border_vals(dirs, omega, eps)
% Determine the border values for the fields.

global DIMS_
dims = DIMS_;

% A simple guesstimation of the length from input to output port.
% Used to guess the initial phase difference between input and output port.
% If the structure is a uniform waveguide, this method should be exact.
if ((dirs{1}(1) == 'x') & (dirs{2}(1) == 'x'))
    len = dims(1);
elseif ((dirs{1}(1) == 'y') & (dirs{2}(1) == 'y'))
    len = dims(2);
else
    len = 0.5 * (dims(1) + dims(2));
end


    %
    %  Find the modes at the appropriate boundaries.
    %

in_out = {'in', 'out'};
for k = 1 : 2
    % Find mode.
    mode = mode_solve(mode_cutout(eps, dirs{k}), omega, dirs{k});

    % Split phase difference between input and output ports.
    phase = -mode.beta * (len-1);
    if strcmp(in_out{k}, 'in')
        phase = phase/2;
    else
        phase = -phase/2;
    end

    % Insert single mode into empty 2D array.
    [Ex{k}, Ey{k}, Hz{k}] = mode_insert(mode, dirs{k}, in_out{k}, phase);
end


    %
    % Combine both input and output modes together.
    %

Ex = Ex{1} + Ex{2};
Ey = Ey{1} + Ey{2};
Hz = Hz{1} + Hz{2};
