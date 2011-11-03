function [delta, M, W, h] = lbfgs_update(x, g, n_max, h)
% s is x_cur - x_old
% y is g_cur - g_old
% If h is empty, then will restart.
% Note: This function can be made much more efficient, since there is a lot of 
% redundant computation going on here.


    %
    % Check if we need to restart.
    %

if isempty(h) 
    % For restart, we simply guess a scaling value of 1.
    delta = 1;
    W = zeros(length(x),0);
    M = zeros(0);
    h.S = [];
    h.Y = [];
    h.x_prev = x;
    h.g_prev = g;
    return
end


    %
    % Add current (s, y) to h, the history of previous iterations.
    %

s = x - h.x_prev;
y = g - h.g_prev;

h.x_prev = x;
h.g_prev = g;

if size(h.S, 2) < n_max % History not full.
    n = size(h.S, 2) + 1;
else % History full, delete oldest entry.
    n = n_max;
    h.S = h.S(:, 2:n_max);
    h.Y = h.Y(:, 2:n_max);
end

% Check curvature condition.
if ((s' * y) <= 0)
    error('Curvature condition broken (sk^T yk = %e)!', (s' * y));
    % TODO: Replace y (damped update) when condition is broken.
end

% Insert new values of s and y into the history.
h.S(:,n) = s;
h.Y(:,n) = y;


    %
    % Form delta, W and M.
    %

delta = dot(y, y) / dot(s, y); % Scaling factor.

% Matrices used to form M.
D = diag(dot(h.S, h.Y, 1));
L = h.S' * h.Y;
L = L - triu(L); % Only keep elements below the diagonal.

% Form M and W.
M = inv([delta*h.S'*h.S, L; L', -D]);
W = [delta*h.S, h.Y];
