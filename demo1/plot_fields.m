function plot_fields(dims, varargin)

N = length(varargin);

if (N > 3)
    xx = ceil(N/2);
    yy = 2;
else
    xx = N;
    yy = 1;
end

for k = 1 : N
    v = varargin{k};
    subplot(yy, xx, k);
    imagesc(real(reshape(v{2}, dims))', 1 * max(abs(v{2}(:)) + eps) * [-1 1]);
    title(v{1});
    axis equal tight;
    set(gca, 'Ydir', 'normal');
end
colormap('jet');
drawnow

