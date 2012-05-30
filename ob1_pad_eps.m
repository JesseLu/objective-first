function [eps] = ob1_pad_eps(eps, pad)

% Expand eps to the full simulation size.
eps = cat(1, repmat(eps(1,:), pad(1), 1), eps, repmat(eps(end,:), pad(2), 1));
eps = cat(2, repmat(eps(:,1), 1, pad(3)), eps, repmat(eps(:,end), 1, pad(4)));
