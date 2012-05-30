function ob1_plot_images(varargin)
% Plot the images for the user to see.
N = length(varargin);

f = factor(N);
f = f(end:-1:1);
n = length(f);

xx = prod(f(1:2:n));
yy = prod(f(2:2:n));

for k = 1 : N
    subplot(xx, yy, k); 
    [im, map] = imread(varargin{k}{1}); 
    if ~isempty(map) % Try to convert to RGB values if needed (mapped data).
        try
            im = idx2rgb(im, map ); 
        catch 
            colormap('jet');
        end
    end
    image(im);
    title(varargin{k}{2});
    axis equal tight; 
end
