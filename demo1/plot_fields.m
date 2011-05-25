function plot_fields(dims, varargin)

N = length(varargin);

f = factor(N);
f = f(end:-1:1);
n = length(f);

xx = prod(f(1:2:n));
yy = prod(f(2:2:n));

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

