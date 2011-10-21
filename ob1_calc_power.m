function [P_out] = ob1_calc_power(spec, Ex, Ey, Hz, t_pml, pad)

dims = size(Ex);
    %
    % Calculate power output to desired mode.
    %
        
% Location for power calculation.
out_pos = min([round(dims(1)-pad(2)/2), dims(1)-t_pml-1]); 

% Project y onto x.
proj = @(x, y) (dot(y(:), x(:)) / norm(x(:))^2) * x(:);

% Calculate the power in the desired output mode.
calcP = @(loc) 0.5 * real(...
                dot(proj(spec.out.Hz, Hz(loc,pad(3)+1:end-pad(4))), ...
                    proj(spec.out.Ey * exp(0.0 * i * spec.out.beta), ...
                        Ey(loc,pad(3)+1:end-pad(4)))));

out_pos = round(dims(1) - pad(2)) : dims(1) - t_pml - 1;
for k = 1 : length(out_pos)
    P_out(k) = calcP(out_pos(k));
end
plot(P_out, '.-')
P_out
pause
P_out = mean(P_out)

% Calculate power leaving a box.
Pbox = @(x,y) dot(Hz(x,y), Ey(x,y));
box_pad = t_pml + 5;
box = [box_pad, dims(1)-box_pad, box_pad, dims(2)-box_pad];
% bottom = 0.5 * real(Pbox(box(1):box(2),box(3)))
% top = 0.5 * real(Pbox(box(1):box(2),box(4)))
% left = 0.5 * real(Pbox(box(1),box(3):box(4)))
right = 0.5 * real(Pbox(box(2),box(3):box(4)))
% Pbox_total = bottom + top + left + right
  
