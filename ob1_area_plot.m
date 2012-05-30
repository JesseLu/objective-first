function ob1_area_plot(z, filename, varargin)
% Normalized area plot.
h = area(1:length(z), z./max(abs(z)));
% set(h, 'FaceColor', [255 194 0]./256); % Tangerine.

if isempty(varargin)
    axis([1 length(z) -3 3]);
elseif strcmp(varargin{1}, 'pos') % Only positive values to be plotted.
    axis([1 length(z) -2 3]);
else
    error('Invalid option.');
end
    
set(gca, 'ytick', []); % No ticks wanted.
print(gcf, '-dpng', '-r150', [filename]); % Save image.
[im] = imread([filename, '.png']); % Reload image.
im = my_add_border(im(300:569,160:1059,:), 0); % Crop and add border.
imwrite(uint8(im), [filename, '.png'], 'png'); % Save.

function [A1] = my_add_border(A0, val)
% Add a one pixel border around image.
A1 = val * ones(size(A0));
A1(2:end-1,2:end-1,:) = A0(2:end-1,2:end-1,:);


