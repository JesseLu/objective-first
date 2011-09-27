function ob1_plot(x, p, dims, option, varargin)

N = prod(dims);

% Ex = reshape(x(1:N), dims);
% Ey = reshape(x(N+1:2*N), dims);
% Hz = reshape(x(2*N+1:3*N), dims);
Hz = reshape(x, dims);

e = reshape(1./p, dims);
% ex = reshape(epsilon(1:N), dims);
% ey = reshape(epsilon(N+1:2*N), dims);

switch option
    case 'quick'
        plot_fields(dims, {'\epsilon', e}, {'Im(Hz)', imag(Hz)}, {'|Hz|', abs(Hz)});
        
    case 'full'
        figure(1);
        plot_fields(dims, {'\epsilon', e});

        figure(2);
        plot_fields(dims, ...
            {'Re(Hz)', real(Hz)}, {'Im(Hz)', imag(Hz)}, {'|Hz|', abs(Hz)});
end

drawnow;

if ~isempty(varargin)
    % Make a picture
    saveas(gcf, ['temp/plot', num2str(1e4+varargin{1}), '.png']);
end
