function waves()
    d = 100;
    omega = .15;
    cutoff = 1 * omega;
    width = 5;

    dims = [200 200];
    y_pos = [1 : dims(2)] - dims(2)/2;
    gauss = @(pos, sigma) exp(-(y_pos - pos).^2 ./ (2*sigma));
    h = truncate_freqs(gauss(0, width), cutoff);

    source = zeros(dims);
    source(170-d,:) = backtrack(omega, h, d, cutoff);
    % source(70-(d-1),:) = backtrack(omega, h, d-1, omega);

    [Ex, Ey, Hz0] = ob1_fdfd_adv(omega, {ones(dims), ones(dims)}, ...
                                ones(dims), source, 'per', 20);
    figure(1)
    ob1_plot(dims, {'|Hz|', abs(Hz0)}, {'real(Hz)', real(Hz0)});
%     z = backtrack(omega,h,d);
%     H = fft(h);
    figure(2);
    subplot 311; plot([abs(Hz0(170,:)); abs(h)]', '.-');
    % imagesc(abs(Hz0)', [0 3]);
%     subplot 312; plot([abs(z); real(z); imag(z)]', '.-');
%     subplot 313; plot([abs(H); real(H); imag(H)]', '.-');

function [source] = backtrack(omega, h, d, kx_cutoff)
    n = length(h);
    kx = fftshift(2*pi*([0:n-1] - n/2) / n);
    ky = sqrt((omega)^2 - kx.^2); % Get ky vectors

    H = fft(h); % Get the plane wave components.
    H = H .* exp(-i * ky * d); % Backtrack a distance d.
    H = H .* (abs(kx) <= kx_cutoff);
    source = ifft(H);


function [z] = truncate_freqs(z, kx_cutoff)
    n = length(z);
    kx = fftshift(2*pi*([0:n-1] - n/2) / n);
    z = ifft((abs(kx) <= kx_cutoff) .* fft(z));
    
