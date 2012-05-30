function ob1_imagesc(z, map, lims, filename, varargin)
% Write out a mapped image.
z = (((z)-lims(1)) / diff(lims) * 63) + 1;
z = 1 * (z < 1) + 64 * (z > 64) + z .* ((z >= 1) & (z <= 64));
if ~isempty(varargin) % Draw a vertical dotted line.
    xline = varargin{1};
    z(xline, 1:2:end) = 1;
end
imwrite(z', map, [filename, '.png']);
imwrite([64:-1:1]', map, [filename, '_cbar.png']); % Colorbar.


