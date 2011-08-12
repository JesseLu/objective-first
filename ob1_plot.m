function ob1_plot(x, epsilon, dims, option)

N = prod(dims);

% Ex = reshape(x(1:N), dims);
% Ey = reshape(x(N+1:2*N), dims);
% Hz = reshape(x(2*N+1:3*N), dims);
Hz = reshape(x, dims);

e = reshape(epsilon, dims);
% ex = reshape(epsilon(1:N), dims);
% ey = reshape(epsilon(N+1:2*N), dims);

switch option
    case 'quick'
        plot_fields(dims, {'\epsilon', e}, {'|Hz|', abs(Hz)});
        
    case 'full'
        figure(1);
        plot_fields(dims, {'\epsilon', e});

        figure(2);
        plot_fields(dims, ...
            {'Re(Hz)', real(Hz)}, {'Im(Hz)', imag(Hz)}, {'|Hz|', abs(Hz)});
end

drawnow;
