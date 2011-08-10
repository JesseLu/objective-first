function ob1_plot(x, epsilon, dims, option)

N = prod(dims);

Ex = reshape(x(1:N), dims);
Ey = reshape(x(N+1:2*N), dims);
Hz = reshape(x(2*N+1:3*N), dims);

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
            {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, ...
            {'Re(Hz)', real(Hz)}, {'Im(Ex)', imag(Ex)}, ...
            {'Im(Ey)', imag(Ey)}, {'Im(Hz)', imag(Hz)}, ...
            {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});
end

drawnow;
