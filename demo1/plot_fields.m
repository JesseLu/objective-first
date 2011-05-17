function plot_fields(dims, varargin)

N = length(varargin);
for k = 1 : N
    v = varargin{k};
    subplot(1, N, k);
    imagesc(real(reshape(v{2}, dims))', mean(abs(v{2}(:)) + eps) * [-1 1]);
    title(v{1});
    axis equal tight;
    set(gca, 'Ydir', 'normal');
end
drawnow

