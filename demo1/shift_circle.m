function [A] = circle_shift(dims, shift)
% A = CIRCLE_SHIFT(DIMS, SHIFT)
%
% Description
%     Create a matrix that will shift the values of a vector, where the vector
%     represents a 2D grid of values. Assumes periodic/circular boundary 
%     conditions.
% 
% Inputs
%     DIMS: 2-element vector of positive integers.
%         DIMS = [XX YY], where XX and YY are the number of grid points in the
%         horizontal and vertical directions, respectively, of a grid of values
%         in vector form.
% 
%     SHIFT: 2-element vector of integers.
%         SHIFT consists of [SX SY], where SX is the horizontal shift and SY 
%         is the vertical shift. Positive values shift downwards and to the 
%         right, while negative values shift upwards and to the left.
% 
% Outputs
%     A: Sparse matrix.
%         Sparse matrix that will perform the shift.


N = prod(dims); % Final matrix will be NxN.


    % 
    % Matrix with place-holder values.
    %

A = reshape(1:N, dims);


    %
    % Shift in the horizontal (x) direction.
    %

if (shift(1) > 0)
    A = [A(end-shift(1)+1:end,:); A(1:end-shift(1),:)];
elseif (shift(1) < 0)
    A = [A(-shift(1)+1:end,:); A(1:-shift(1),:)];
end


    %
    % Shift in the vertical (y) direction.
    %

if (shift(2) > 0)
    A = [A(:,end-shift(2)+1:end), A(:,1:end-shift(2))];
elseif (shift(2) < 0)
    A = [A(:,-shift(2)+1:end), A(:,1:-shift(2))];
end


    %
    % Now form the actual, sparse shifting matrix.
    %

A = sparse(1:N, A(:), ones(1, N), N, N);
